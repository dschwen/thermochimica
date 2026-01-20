#include "thermochimica/solver/GEMSolver.hpp"
#include <Eigen/LU>
#include <cmath>
#include <iostream>

namespace Thermochimica {

int GEMNewton::compute(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;
    auto& gem = *ctx.gem;

    // Determine system size
    int nVar = thermo.nElements + thermo.nConPhases + thermo.nSolnPhases;

    // Construct Hessian and RHS
    Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(nVar, nVar);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(nVar);

    constructHessian(ctx, hessian, rhs);

    // Solve using Eigen's PartialPivLU (equivalent to LAPACK DGESV)
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(hessian);

    // Check for singularity
    if (std::abs(lu.determinant()) < 1e-300) {
        handleSingular(ctx, hessian, rhs);
        lu.compute(hessian);
    }

    // Solve the system
    gem.dUpdateVar.head(nVar) = lu.solve(rhs);

    // Check for NaN in solution
    for (int i = 0; i < nVar; ++i) {
        if (std::isnan(gem.dUpdateVar(i)) || std::isinf(gem.dUpdateVar(i))) {
            gem.lRevertSystem = true;
            return 1;
        }
    }

    return 0;
}

void GEMNewton::constructHessian(ThermoContext& ctx,
                                  Eigen::MatrixXd& hessian,
                                  Eigen::VectorXd& rhs) {
    auto& thermo = *ctx.thermo;

    int nElements = thermo.nElements;
    int nConPhases = thermo.nConPhases;
    int nSolnPhases = thermo.nSolnPhases;

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

    // Construct RHS (constraint vector)
    // Element mass balance constraints
    for (int j = 0; j < nElements; ++j) {
        double sum = thermo.dMolesElement(j);

        // Subtract contribution from solution species
        for (int iPhase = 0; iPhase < nSolnPhases; ++iPhase) {
            int phaseIdx = -thermo.iAssemblage(nElements - nSolnPhases + iPhase) - 1;
            if (phaseIdx < 0) continue;

            int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
            int iLast = thermo.nSpeciesPhase(phaseIdx + 1);

            for (int iSpecies = iFirst; iSpecies < iLast; ++iSpecies) {
                if (thermo.dMolFraction(iSpecies) > 0) {
                    sum += (thermo.dChemicalPotential(iSpecies) - 1.0) *
                           thermo.dMolesSpecies(iSpecies) *
                           thermo.dStoichSpecies(iSpecies, j) /
                           thermo.iParticlesPerMole(iSpecies);
                }
            }
        }

        rhs(j) = sum;
    }

    // Solution phase Gibbs energy constraints
    for (int iPhase = 0; iPhase < nSolnPhases; ++iPhase) {
        int row = nElements + iPhase;
        int phaseIdx = -thermo.iAssemblage(nElements - nSolnPhases + iPhase) - 1;
        if (phaseIdx < 0) continue;

        rhs(row) = thermo.dGibbsSolnPhase(phaseIdx);
    }

    // Pure condensed phase constraints
    for (int iPhase = 0; iPhase < nConPhases; ++iPhase) {
        int row = nElements + nSolnPhases + iPhase;
        int speciesIdx = thermo.iAssemblage(iPhase) - 1;
        if (speciesIdx < 0 || speciesIdx >= thermo.nSpecies) continue;

        rhs(row) = thermo.dStdGibbsEnergy(speciesIdx);
    }
}

void GEMNewton::handleSingular(ThermoContext& ctx,
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
