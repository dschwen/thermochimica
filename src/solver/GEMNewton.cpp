#include "thermochimica/solver/GEMSolver.hpp"
#include "thermochimica/context/PhaseConstraints.hpp"
#include "thermochimica/util/Constants.hpp"
#include <Eigen/LU>
#include <cmath>
#include <iostream>

namespace Thermochimica {

int GEMNewton::compute(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;
    auto& gem = *ctx.gem;

    // Build list of active elements (with non-zero input mass)
    std::vector<int> activeElements;
    for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
        if (thermo.dMolesElement(j) > 0.0) {
            activeElements.push_back(j);
        }
    }
    int nActiveEl = static_cast<int>(activeElements.size());

    // Determine reduced system size (active elements + phases)
    int nVar = nActiveEl + thermo.nConPhases + thermo.nSolnPhases;

    // Construct Hessian and RHS for reduced system
    Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(nVar, nVar);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(nVar);

    constructHessianReduced(ctx, hessian, rhs, activeElements);

    // Add phase constraint penalty contributions to Newton system
    if (ctx.phaseConstraints->hasActiveConstraints()) {
        ConstrainedGEM::addConstraintGradients(ctx, rhs, activeElements);
        ConstrainedGEM::addConstraintHessian(ctx, hessian, activeElements);
    }

    static int debugNewtonCount = 0;

    // Debug: print Hessian and RHS before solving
    if (debugNewtonCount < 2) {
        std::cerr << "[GEMNewton] Reduced Hessian (" << nVar << "x" << nVar << ") with " << nActiveEl << " active elements:\n";
        for (int i = 0; i < nVar; ++i) {
            std::cerr << "  ";
            for (int j = 0; j < nVar; ++j) {
                std::cerr << hessian(i, j) << " ";
            }
            std::cerr << "\n";
        }
        std::cerr << "[GEMNewton] RHS: ";
        for (int i = 0; i < nVar; ++i) {
            std::cerr << rhs(i) << " ";
        }
        std::cerr << "\n";
    }

    // Solve using Eigen's PartialPivLU (equivalent to LAPACK DGESV)
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(hessian);

    // Check for singularity
    if (std::abs(lu.determinant()) < 1e-20) {
        handleSingularReduced(hessian, rhs);
        lu.compute(hessian);
    }

    // Solve the reduced system
    Eigen::VectorXd updateReduced = lu.solve(rhs);

    // Map back to full update vector
    gem.dUpdateVar.setZero();
    for (int i = 0; i < nActiveEl; ++i) {
        gem.dUpdateVar(activeElements[i]) = updateReduced(i);
    }
    // Solution phase updates (map from reduced to full indexing)
    for (int i = 0; i < thermo.nSolnPhases; ++i) {
        gem.dUpdateVar(thermo.nElements + i) = updateReduced(nActiveEl + i);
    }
    // Condensed phase updates
    for (int i = 0; i < thermo.nConPhases; ++i) {
        gem.dUpdateVar(thermo.nElements + thermo.nSolnPhases + i) = updateReduced(nActiveEl + thermo.nSolnPhases + i);
    }

    // Check for NaN in solution
    int nVarFull = thermo.nElements + thermo.nConPhases + thermo.nSolnPhases;
    for (int i = 0; i < nVarFull; ++i) {
        if (std::isnan(gem.dUpdateVar(i)) || std::isinf(gem.dUpdateVar(i))) {
            gem.lRevertSystem = true;
            return 1;
        }
        if (debugNewtonCount < 2) {
            std::cerr << "[GEMNewton] dUpdateVar(" << i << ")=" << gem.dUpdateVar(i) << "\n";
        }
    }
    if (debugNewtonCount < 2) {
        std::cerr << "[GEMNewton] det=" << lu.determinant() << " nVar=" << nVar << "\n";
    }
    ++debugNewtonCount;

    return 0;
}

void GEMNewton::constructHessianReduced(ThermoContext& ctx,
                                         Eigen::MatrixXd& hessian,
                                         Eigen::VectorXd& rhs,
                                         const std::vector<int>& activeElements) {
    auto& thermo = *ctx.thermo;
    auto& io = *ctx.io;

    int nActiveEl = static_cast<int>(activeElements.size());
    int nConPhases = thermo.nConPhases;
    int nSolnPhases = thermo.nSolnPhases;
    double R = Constants::kIdealGasConstant;
    double T = io.dTemperature;

    // Build Hessian matrix blocks (using reduced element indexing)
    // Block 1: Element-element block (from solution phase species)
    for (int iPhase = 0; iPhase < nSolnPhases; ++iPhase) {
        int phaseIdx = -thermo.iAssemblage(thermo.nElements - nSolnPhases + iPhase) - 1;
        if (phaseIdx < 0) continue;

        int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
        int iLast = thermo.nSpeciesPhase(phaseIdx + 1);

        for (int iSpecies = iFirst; iSpecies < iLast; ++iSpecies) {
            if (thermo.dMolFraction(iSpecies) <= 0) continue;

            double molesSpec = thermo.dMolesSpecies(iSpecies);
            int particles = thermo.iParticlesPerMole(iSpecies);
            double factor = molesSpec / (particles * particles);

            for (int jRed = 0; jRed < nActiveEl; ++jRed) {
                int jFull = activeElements[jRed];
                for (int iRed = jRed; iRed < nActiveEl; ++iRed) {
                    int iFull = activeElements[iRed];
                    double contrib = thermo.dStoichSpecies(iSpecies, iFull) *
                                    thermo.dStoichSpecies(iSpecies, jFull) * factor;
                    hessian(iRed, jRed) += contrib;
                    if (iRed != jRed) {
                        hessian(jRed, iRed) += contrib;  // Symmetry
                    }
                }
            }
        }
    }

    // Block 2: Element-Solution phase block
    for (int iPhase = 0; iPhase < nSolnPhases; ++iPhase) {
        int col = nActiveEl + iPhase;
        int phaseIdx = -thermo.iAssemblage(thermo.nElements - nSolnPhases + iPhase) - 1;
        if (phaseIdx < 0) continue;

        int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
        int iLast = thermo.nSpeciesPhase(phaseIdx + 1);

        for (int jRed = 0; jRed < nActiveEl; ++jRed) {
            int jFull = activeElements[jRed];
            double sum = 0.0;
            for (int iSpecies = iFirst; iSpecies < iLast; ++iSpecies) {
                if (thermo.dMolFraction(iSpecies) > 0) {
                    sum += thermo.dStoichSpecies(iSpecies, jFull) *
                           thermo.dMolFraction(iSpecies) /
                           thermo.iParticlesPerMole(iSpecies);
                }
            }
            hessian(jRed, col) = sum;
            hessian(col, jRed) = sum;  // Symmetry
        }
    }

    // Block 3: Element-Condensed phase block
    for (int iPhase = 0; iPhase < nConPhases; ++iPhase) {
        int col = nActiveEl + nSolnPhases + iPhase;
        int speciesIdx = thermo.iAssemblage(iPhase) - 1;
        if (speciesIdx < 0 || speciesIdx >= thermo.nSpecies) continue;

        for (int jRed = 0; jRed < nActiveEl; ++jRed) {
            int jFull = activeElements[jRed];
            double val = thermo.dStoichSpecies(speciesIdx, jFull) /
                        thermo.iParticlesPerMole(speciesIdx);
            hessian(jRed, col) = val;
            hessian(col, jRed) = val;
        }
    }

    // Construct RHS as residuals (negative gradient for Newton)
    // Newton solves: H * Δx = -r, where r is the residual

    // Element mass balance residuals (for active elements only)
    // Residual = target - current = b_j - ∑_i (n_i * a_ij / p_i)
    for (int jRed = 0; jRed < nActiveEl; ++jRed) {
        int jFull = activeElements[jRed];
        double target = thermo.dMolesElement(jFull);
        double current = 0.0;

        // Sum contribution from solution phase species
        for (int iPhase = 0; iPhase < nSolnPhases; ++iPhase) {
            int phaseIdx = -thermo.iAssemblage(thermo.nElements - nSolnPhases + iPhase) - 1;
            if (phaseIdx < 0) continue;

            int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
            int iLast = thermo.nSpeciesPhase(phaseIdx + 1);

            for (int iSpecies = iFirst; iSpecies < iLast; ++iSpecies) {
                if (thermo.dMolFraction(iSpecies) > 0) {
                    current += thermo.dMolesSpecies(iSpecies) *
                               thermo.dStoichSpecies(iSpecies, jFull) /
                               thermo.iParticlesPerMole(iSpecies);
                }
            }
        }

        // Sum contribution from pure condensed phase species
        for (int iPhase = 0; iPhase < nConPhases; ++iPhase) {
            int speciesIdx = thermo.iAssemblage(iPhase) - 1;
            if (speciesIdx < 0 || speciesIdx >= thermo.nSpecies) continue;

            current += thermo.dMolesPhase(iPhase) *
                       thermo.dStoichSpecies(speciesIdx, jFull) /
                       thermo.iParticlesPerMole(speciesIdx);
        }

        // RHS = negative residual for Newton method
        rhs(jRed) = -(target - current);
    }

    // Solution phase Gibbs energy residuals
    // At equilibrium: dGibbsSolnPhase = 0 (normalized driving force)
    // Residual = 0 - current = -dGibbsSolnPhase
    for (int iPhase = 0; iPhase < nSolnPhases; ++iPhase) {
        int row = nActiveEl + iPhase;
        int phaseIdx = -thermo.iAssemblage(thermo.nElements - nSolnPhases + iPhase) - 1;
        if (phaseIdx < 0) continue;

        // RHS = negative residual = -(-dGibbsSolnPhase) = dGibbsSolnPhase
        rhs(row) = thermo.dGibbsSolnPhase(phaseIdx);
    }

    // Pure condensed phase equilibrium residuals
    // At equilibrium: G°/RT = ∑_j λ_j * a_ij / p_i
    // Residual = G°/RT - ∑_j λ_j * a_ij / p_i
    for (int iPhase = 0; iPhase < nConPhases; ++iPhase) {
        int row = nActiveEl + nSolnPhases + iPhase;
        int speciesIdx = thermo.iAssemblage(iPhase) - 1;
        if (speciesIdx < 0 || speciesIdx >= thermo.nSpecies) continue;

        double gRT = thermo.dStdGibbsEnergy(speciesIdx) / (R * T);

        // Compute current ∑_j λ_j * a_ij / p_i
        double lambdaSum = 0.0;
        for (int jRed = 0; jRed < nActiveEl; ++jRed) {
            int jFull = activeElements[jRed];
            lambdaSum += thermo.dElementPotential(jFull) *
                         thermo.dStoichSpecies(speciesIdx, jFull);
        }
        lambdaSum /= thermo.iParticlesPerMole(speciesIdx);

        // RHS = negative residual for Newton method
        rhs(row) = -(gRT - lambdaSum);
    }
}

void GEMNewton::handleSingularReduced(Eigen::MatrixXd& hessian, Eigen::VectorXd& rhs) {
    // For singular matrices, add small diagonal perturbation
    int n = hessian.rows();
    for (int i = 0; i < n; ++i) {
        if (std::abs(hessian(i, i)) < 1e-20) {
            hessian(i, i) = 1.0;
            rhs(i) = 0.0;
        }
    }
}

void GEMNewton::constructHessian(ThermoContext& ctx,
                                  Eigen::MatrixXd& hessian,
                                  Eigen::VectorXd& rhs) {
    auto& thermo = *ctx.thermo;
    auto& io = *ctx.io;

    int nElements = thermo.nElements;
    int nConPhases = thermo.nConPhases;
    int nSolnPhases = thermo.nSolnPhases;
    double R = Constants::kIdealGasConstant;
    double T = io.dTemperature;

    // Build Hessian matrix blocks
    // Block 1: Element-element block (from solution phase species)
    for (int iPhase = 0; iPhase < nSolnPhases; ++iPhase) {
        int phaseIdx = -thermo.iAssemblage(nElements - nSolnPhases + iPhase) - 1;
        if (phaseIdx < 0) continue;

        int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
        int iLast = thermo.nSpeciesPhase(phaseIdx + 1);

        for (int iSpecies = iFirst; iSpecies < iLast; ++iSpecies) {
            if (thermo.dMolFraction(iSpecies) <= 0) continue;

            double molesSpec = thermo.dMolesSpecies(iSpecies);
            int particles = thermo.iParticlesPerMole(iSpecies);
            double factor = molesSpec / (particles * particles);

            for (int j = 0; j < nElements; ++j) {
                for (int i = j; i < nElements; ++i) {
                    double contrib = thermo.dStoichSpecies(iSpecies, i) *
                                    thermo.dStoichSpecies(iSpecies, j) * factor;
                    hessian(i, j) += contrib;
                    if (i != j) {
                        hessian(j, i) += contrib;  // Symmetry
                    }
                }
            }
        }
    }

    // Block 2: Element-Solution phase block
    for (int iPhase = 0; iPhase < nSolnPhases; ++iPhase) {
        int col = nElements + iPhase;
        int phaseIdx = -thermo.iAssemblage(nElements - nSolnPhases + iPhase) - 1;
        if (phaseIdx < 0) continue;

        int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
        int iLast = thermo.nSpeciesPhase(phaseIdx + 1);

        for (int j = 0; j < nElements; ++j) {
            double sum = 0.0;
            for (int iSpecies = iFirst; iSpecies < iLast; ++iSpecies) {
                if (thermo.dMolFraction(iSpecies) > 0) {
                    sum += thermo.dStoichSpecies(iSpecies, j) *
                           thermo.dMolFraction(iSpecies) /
                           thermo.iParticlesPerMole(iSpecies);
                }
            }
            hessian(j, col) = sum;
            hessian(col, j) = sum;  // Symmetry
        }
    }

    // Block 3: Element-Condensed phase block
    for (int iPhase = 0; iPhase < nConPhases; ++iPhase) {
        int col = nElements + nSolnPhases + iPhase;
        int speciesIdx = thermo.iAssemblage(iPhase) - 1;
        if (speciesIdx < 0 || speciesIdx >= thermo.nSpecies) continue;

        for (int j = 0; j < nElements; ++j) {
            double val = thermo.dStoichSpecies(speciesIdx, j) /
                        thermo.iParticlesPerMole(speciesIdx);
            hessian(j, col) = val;
            hessian(col, j) = val;
        }
    }

    // Construct RHS as residuals (negative gradient for Newton)
    // Newton solves: H * Δx = -r, where r is the residual

    // Element mass balance residuals
    // Residual = target - current = b_j - ∑_i (n_i * a_ij / p_i)
    for (int j = 0; j < nElements; ++j) {
        double target = thermo.dMolesElement(j);
        double current = 0.0;

        // Sum contribution from solution phase species
        for (int iPhase = 0; iPhase < nSolnPhases; ++iPhase) {
            int phaseIdx = -thermo.iAssemblage(nElements - nSolnPhases + iPhase) - 1;
            if (phaseIdx < 0) continue;

            int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
            int iLast = thermo.nSpeciesPhase(phaseIdx + 1);

            for (int iSpecies = iFirst; iSpecies < iLast; ++iSpecies) {
                if (thermo.dMolFraction(iSpecies) > 0) {
                    current += thermo.dMolesSpecies(iSpecies) *
                               thermo.dStoichSpecies(iSpecies, j) /
                               thermo.iParticlesPerMole(iSpecies);
                }
            }
        }

        // Sum contribution from pure condensed phase species
        for (int iPhase = 0; iPhase < nConPhases; ++iPhase) {
            int speciesIdx = thermo.iAssemblage(iPhase) - 1;
            if (speciesIdx < 0 || speciesIdx >= thermo.nSpecies) continue;

            current += thermo.dMolesPhase(iPhase) *
                       thermo.dStoichSpecies(speciesIdx, j) /
                       thermo.iParticlesPerMole(speciesIdx);
        }

        // RHS = negative residual for Newton method
        rhs(j) = -(target - current);
    }

    // Solution phase Gibbs energy residuals
    // At equilibrium: dGibbsSolnPhase = 0
    for (int iPhase = 0; iPhase < nSolnPhases; ++iPhase) {
        int row = nElements + iPhase;
        int phaseIdx = -thermo.iAssemblage(nElements - nSolnPhases + iPhase) - 1;
        if (phaseIdx < 0) continue;

        rhs(row) = thermo.dGibbsSolnPhase(phaseIdx);
    }

    // Pure condensed phase equilibrium residuals
    // At equilibrium: G°/RT = ∑_j λ_j * a_ij / p_i
    for (int iPhase = 0; iPhase < nConPhases; ++iPhase) {
        int row = nElements + nSolnPhases + iPhase;
        int speciesIdx = thermo.iAssemblage(iPhase) - 1;
        if (speciesIdx < 0 || speciesIdx >= thermo.nSpecies) continue;

        double gRT = thermo.dStdGibbsEnergy(speciesIdx) / (R * T);

        // Compute current ∑_j λ_j * a_ij / p_i
        double lambdaSum = 0.0;
        for (int j = 0; j < nElements; ++j) {
            lambdaSum += thermo.dElementPotential(j) *
                         thermo.dStoichSpecies(speciesIdx, j);
        }
        lambdaSum /= thermo.iParticlesPerMole(speciesIdx);

        // RHS = negative residual for Newton method
        rhs(row) = -(gRT - lambdaSum);
    }
}

void GEMNewton::handleSingular(ThermoContext& /*ctx*/,
                                Eigen::MatrixXd& hessian,
                                Eigen::VectorXd& rhs) {
    // For singular matrices, add small diagonal perturbation
    int n = hessian.rows();
    for (int i = 0; i < n; ++i) {
        if (std::abs(hessian(i, i)) < 1e-20) {
            hessian(i, i) = 1.0;
            rhs(i) = 0.0;
        }
    }
}

} // namespace Thermochimica
