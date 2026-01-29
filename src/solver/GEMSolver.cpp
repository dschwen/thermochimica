#include "thermochimica/solver/GEMSolver.hpp"
#include "thermochimica/context/PhaseConstraints.hpp"
#include "thermochimica/util/ErrorCodes.hpp"
#include "thermochimica/util/Constants.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace Thermochimica {

// Forward declarations
void computeExcessGibbs(ThermoContext& ctx, int phaseIndex);

/// @brief Compute chemical potentials for all solution phases
static void computeChemicalPotentials(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;
    auto& gem = *ctx.gem;

    // For each stable solution phase
    for (int iPhase = 0; iPhase < thermo.nSolnPhases; ++iPhase) {
        // Get the phase index from assemblage
        int assembIdx = thermo.nElements - thermo.nSolnPhases + iPhase;
        if (assembIdx < 0 || assembIdx >= static_cast<int>(thermo.iAssemblage.size())) {
            continue;
        }
        int phaseIdx = -thermo.iAssemblage(assembIdx) - 1;

        if (phaseIdx < 0 || phaseIdx >= thermo.nSolnPhasesSys) {
            continue;
        }

        int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
        int iLast = thermo.nSpeciesPhase(phaseIdx + 1);

        // Reset chemical potentials and excess Gibbs contributions
        for (int i = iFirst; i < iLast; ++i) {
            // Initialize with standard state (already computed in compThermoData)
            // thermo.dChemicalPotential is already set to G0/RT
            // Reset excess contributions
            gem.dPartialExcessGibbs(i) = 0.0;
        }

        // Compute excess Gibbs energy contributions
        computeExcessGibbs(ctx, phaseIdx);

        // Add ideal mixing and excess contributions to chemical potential
        // Note: dStdGibbsEnergy is in J/mol, need to normalize by RT to get dimensionless mu
        double R = Constants::kIdealGasConstant;
        double T = ctx.io->dTemperature;
        for (int i = iFirst; i < iLast; ++i) {
            double x = thermo.dMolFraction(i);
            double muStd = thermo.dStdGibbsEnergy(i) / (R * T);  // Dimensionless
            if (x > 1e-100) {
                // Chemical potential = G°/RT + ln(x) + excess/RT
                thermo.dChemicalPotential(i) = muStd + std::log(x);
                // Add excess contribution (in J/mol, need to normalize by RT)
                thermo.dChemicalPotential(i) += gem.dPartialExcessGibbs(i) / (R * T);
            } else {
                thermo.dChemicalPotential(i) = muStd + std::log(1e-100);
            }
        }

        // Compute driving force (dGibbsSolnPhase) for this phase
        // DF = ∑_i x_i * (μ_i - μ*_i)
        // where μ*_i = ∑_j λ_j * a_{ij} / p_i
        // At equilibrium, DF = 0
        double drivingForce = 0.0;
        for (int i = iFirst; i < iLast; ++i) {
            double x = thermo.dMolFraction(i);
            if (x > 1e-100) {
                // Compute mu* from element potentials
                double muStar = 0.0;
                for (int j = 0; j < thermo.nElements; ++j) {
                    muStar += thermo.dElementPotential(j) * thermo.dStoichSpecies(i, j);
                }
                int particles = thermo.iParticlesPerMole(i);
                muStar /= static_cast<double>(particles);

                // Driving force contribution
                drivingForce += x * (thermo.dChemicalPotential(i) - muStar);
            }
        }
        thermo.dGibbsSolnPhase(phaseIdx) = drivingForce;
    }
}

/// Run inner GEM iteration loop (standard or with constraints)
/// @return true if converged
static bool runInnerGEMLoop(ThermoContext& ctx) {
    auto& io = *ctx.io;
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;
    auto& pc = *ctx.phaseConstraints;

    bool hasConstraints = pc.hasActiveConstraints();

    // Main inner iteration loop
    for (gem.iterGlobal = 1; gem.iterGlobal <= Constants::kIterGlobalMax; ++gem.iterGlobal) {
        // Compute chemical potentials for solution phases (includes excess Gibbs)
        if (thermo.nSolnPhases > 0) {
            computeChemicalPotentials(ctx);
        }

        // Ensure species moles are computed from phase moles and mole fractions
        // This is needed for the Newton step to build a proper Hessian
        for (int iPhase = 0; iPhase < thermo.nSolnPhases; ++iPhase) {
            int assembIdx = thermo.nElements - thermo.nSolnPhases + iPhase;
            int phaseIdx = -thermo.iAssemblage(assembIdx) - 1;
            if (phaseIdx < 0 || phaseIdx >= thermo.nSolnPhasesSys) continue;

            int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
            int iLast = thermo.nSpeciesPhase(phaseIdx + 1);
            double phaseMoles = thermo.dMolesPhase(assembIdx);

            for (int i = iFirst; i < iLast; ++i) {
                if (thermo.dMolFraction(i) > 0) {
                    thermo.dMolesSpecies(i) = phaseMoles * thermo.dMolFraction(i);
                }
            }
        }

        // Also update species moles for pure condensed phases
        for (int i = 0; i < thermo.nConPhases; ++i) {
            int speciesIdx = thermo.iAssemblage(i) - 1;  // Convert from 1-based
            if (speciesIdx >= 0 && speciesIdx < thermo.nSpecies) {
                thermo.dMolesSpecies(speciesIdx) = thermo.dMolesPhase(i);
            }
        }

        // Update mole fractions for IDEAL solution phases based on element potentials
        // At equilibrium: x_i ∝ exp(μ* - μ_std) where μ* = Σ_j λ_j * a_{ij} / p_i
        // Only apply to IDMX phases - non-ideal phases need proper Newton iteration
        if (thermo.nSolnPhases > 0) {
            double R = Constants::kIdealGasConstant;
            double T = io.dTemperature;

            for (int iPhase = 0; iPhase < thermo.nSolnPhases; ++iPhase) {
                int assembIdx = thermo.nElements - thermo.nSolnPhases + iPhase;
                int phaseIdx = -thermo.iAssemblage(assembIdx) - 1;
                if (phaseIdx < 0 || phaseIdx >= thermo.nSolnPhasesSys) continue;

                // Only update mole fractions for ideal mixing phases
                if (phaseIdx >= static_cast<int>(thermo.iSolnPhaseType.size()) ||
                    thermo.iSolnPhaseType[phaseIdx] != Constants::PhaseType::IDMX) {
                    continue;
                }

                int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
                int iLast = thermo.nSpeciesPhase(phaseIdx + 1);

                // Compute new mole fractions from element potentials
                double sumExp = 0.0;
                std::vector<double> expVal(iLast - iFirst, 0.0);

                for (int i = iFirst; i < iLast; ++i) {
                    // Check feasibility
                    bool feasible = true;
                    for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
                        if (thermo.dStoichSpecies(i, j) > 0.0 && thermo.dMolesElement(j) <= 0.0) {
                            feasible = false;
                            break;
                        }
                    }
                    if (!feasible) continue;

                    double muStar = 0.0;
                    for (int j = 0; j < thermo.nElements; ++j) {
                        muStar += thermo.dElementPotential(j) * thermo.dStoichSpecies(i, j);
                    }
                    muStar /= static_cast<double>(thermo.iParticlesPerMole(i));

                    double muStd = thermo.dStdGibbsEnergy(i) / (R * T);
                    double arg = muStar - muStd;
                    // Clamp to prevent overflow
                    arg = std::max(-100.0, std::min(100.0, arg));
                    expVal[i - iFirst] = std::exp(arg);
                    sumExp += expVal[i - iFirst];
                }

                // Update mole fractions with damping
                if (sumExp > 1e-300) {
                    double damping = 0.5;  // Damping factor to prevent oscillation
                    for (int i = iFirst; i < iLast; ++i) {
                        double xNew = expVal[i - iFirst] / sumExp;
                        double xOld = thermo.dMolFraction(i);
                        thermo.dMolFraction(i) = xOld + damping * (xNew - xOld);
                    }

                    // Normalize to ensure sum = 1
                    double sum = 0.0;
                    for (int i = iFirst; i < iLast; ++i) {
                        sum += thermo.dMolFraction(i);
                    }
                    if (sum > 1e-300) {
                        for (int i = iFirst; i < iLast; ++i) {
                            thermo.dMolFraction(i) /= sum;
                        }
                    }
                }
            }
        }

        // For systems with solution phases, compute Newton direction and perform line search
        // Note: Skip Newton when nSolnPhases=0 because the Hessian is singular without
        // solution phase contributions. Phase assemblage will handle adding solution phases.
        //
        // IMPORTANT: Skip Newton for now as it's causing phase moles to explode.
        // The simplified single-phase approach below handles the common case.
        if (false && thermo.nSolnPhases > 0 && thermo.nConPhases > 0) {
            // Full Newton when we have both solution and condensed phases
            GEMNewton::compute(ctx);
            GEMLineSearch::search(ctx);
        } else if (thermo.nSolnPhases >= 1 && thermo.nConPhases == 0) {
            // Solution phases only, no condensed phases: simpler update
            // First, recalculate phase moles to satisfy mass balance
            for (int iPhase = 0; iPhase < thermo.nSolnPhases; ++iPhase) {
                int assembIdx = thermo.nElements - thermo.nSolnPhases + iPhase;
                int phaseIdx = -thermo.iAssemblage(assembIdx) - 1;
                if (phaseIdx < 0 || phaseIdx >= thermo.nSolnPhasesSys) continue;

                int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
                int iLast = thermo.nSpeciesPhase(phaseIdx + 1);

                // Calculate phase moles from mass balance
                // For each active element j: n_φ * Σ_i x_i * a_ij = b_j
                // Use least squares if overdetermined
                double sumNumer = 0.0;
                double sumDenom = 0.0;
                for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
                    if (thermo.dMolesElement(j) <= 0) continue;

                    double avgStoich = 0.0;
                    for (int i = iFirst; i < iLast; ++i) {
                        if (thermo.dMolFraction(i) > 0) {
                            avgStoich += thermo.dMolFraction(i) * thermo.dStoichSpecies(i, j);
                        }
                    }
                    if (avgStoich > 1e-20) {
                        // n_φ = b_j / avgStoich
                        double n_phi_j = thermo.dMolesElement(j) / avgStoich;
                        sumNumer += n_phi_j;
                        sumDenom += 1.0;
                    }
                }

                if (sumDenom > 0) {
                    thermo.dMolesPhase(assembIdx) = sumNumer / sumDenom;
                }

                // Update species moles
                double phaseMoles = thermo.dMolesPhase(assembIdx);
                for (int i = iFirst; i < iLast; ++i) {
                    if (thermo.dMolFraction(i) > 0) {
                        thermo.dMolesSpecies(i) = phaseMoles * thermo.dMolFraction(i);
                    }
                }
            }

            // Now compute element potentials from chemical equilibrium
            // For a single solution phase, element potentials can be computed directly
            // from chemical potentials of pure end-member species
            for (int iPhase = 0; iPhase < thermo.nSolnPhases; ++iPhase) {
                int assembIdx = thermo.nElements - thermo.nSolnPhases + iPhase;
                int phaseIdx = -thermo.iAssemblage(assembIdx) - 1;
                if (phaseIdx < 0 || phaseIdx >= thermo.nSolnPhasesSys) continue;

                int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
                int iLast = thermo.nSpeciesPhase(phaseIdx + 1);

                // Solve for element potentials from chemical equilibrium
                // At equilibrium: μ_i = Σ_j λ_j * a_{ij} / p_i
                for (int i = iFirst; i < iLast; ++i) {
                    if (thermo.dMolFraction(i) < 1e-20) continue;

                    // Check if this is a pure end member (only one element)
                    int singleElement = -1;
                    int nElements = 0;
                    for (int j = 0; j < thermo.nElements; ++j) {
                        if (thermo.dStoichSpecies(i, j) > 0) {
                            singleElement = j;
                            nElements++;
                        }
                    }

                    // For pure end member species, λ_j = μ_i * p_i / a_{ij}
                    if (nElements == 1 && singleElement >= 0) {
                        double stoich = thermo.dStoichSpecies(i, singleElement);
                        double particles = thermo.iParticlesPerMole(i);
                        thermo.dElementPotential(singleElement) =
                            thermo.dChemicalPotential(i) * particles / stoich;
                    }
                }
            }
        }

        // Check phase assemblage
        PhaseAssemblage::check(ctx);

        // If constraints are active, compute current phase fractions
        if (hasConstraints) {
            ConstrainedGEM::computePhaseElementFractions(ctx);
        }

        // Check convergence
        gem.lConverged = ConvergenceChecker::check(ctx);

        // Exit if error or converged
        if (io.INFOThermo != 0 || gem.lConverged) {
            break;
        }
    }

    return gem.lConverged;
}

int GEMSolver::solve(ThermoContext& ctx) {
    auto& io = *ctx.io;
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;
    auto& pc = *ctx.phaseConstraints;

    // Initialize the GEM solver
    init(ctx);

    // Check if system with only pure condensed phases is already converged
    if (thermo.nSolnPhases == 0 && io.INFOThermo == 0 && !io.lReinitLoaded) {
        checkSysOnlyPureConPhases(ctx);
    }

    // Check if phase constraints are active
    if (!pc.hasActiveConstraints()) {
        // Standard unconstrained GEM solve
        bool converged = runInnerGEMLoop(ctx);

        if (!converged && io.INFOThermo == 0) {
            io.INFOThermo = ErrorCode::kGEMSolverDidNotConverge;
        }

        return io.INFOThermo;
    }

    // =========================================================================
    // Augmented Lagrangian outer loop for constrained GEM
    // =========================================================================

    // Reset constraint state
    pc.reset();

    for (pc.currentOuterIteration = 0;
         pc.currentOuterIteration < pc.maxOuterIterations;
         ++pc.currentOuterIteration) {

        // Reset inner solver state for new outer iteration
        gem.reset();
        gem.lConverged = false;

        // Run inner GEM loop with current Lagrange multipliers and penalty
        bool innerConverged = runInnerGEMLoop(ctx);

        if (io.INFOThermo != 0) {
            // Error in inner loop
            return io.INFOThermo;
        }

        if (!innerConverged) {
            // Inner loop didn't converge - try continuing with current result
            // The constraints might still improve
        }

        // Compute current phase element fractions
        ConstrainedGEM::computePhaseElementFractions(ctx);

        // Check if constraints are satisfied
        if (pc.areConstraintsSatisfied()) {
            // All constraints satisfied - we're done
            gem.lConverged = true;
            break;
        }

        // Update Lagrange multipliers: λ += ρ * (f - f_target)
        ConstrainedGEM::updateLagrangeMultipliers(ctx);

        // Increase penalty parameter for next outer iteration
        pc.penaltyParameter *= pc.penaltyGrowthRate;
    }

    // Check final convergence status
    if (!gem.lConverged && io.INFOThermo == 0) {
        // Check if constraints are at least approximately satisfied
        if (!pc.areConstraintsSatisfied()) {
            io.INFOThermo = ErrorCode::kGEMSolverDidNotConverge;
        }
    }

    return io.INFOThermo;
}

void GEMSolver::init(ThermoContext& ctx) {
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;

    // Reset GEM solver state
    gem.reset();

    // Initialize iteration history
    if (gem.iterHistory.rows() < thermo.nElements ||
        gem.iterHistory.cols() < Constants::kIterGlobalMax) {
        gem.iterHistory.resize(thermo.nElements, Constants::kIterGlobalMax);
    }
    gem.iterHistory.setZero();

    // Store initial assemblage in history
    for (int i = 0; i < thermo.nElements; ++i) {
        gem.iterHistory(i, 0) = thermo.iAssemblage(i);
    }

    // Initialize moles phase last
    if (gem.dMolesPhaseLast.size() != thermo.dMolesPhase.size()) {
        gem.dMolesPhaseLast.resize(thermo.dMolesPhase.size());
    }
    gem.dMolesPhaseLast = thermo.dMolesPhase;

    // Initialize solution phase flags
    gem.lSolnPhases.resize(thermo.nSolnPhasesSys, false);
    gem.lMiscibility.resize(thermo.nSolnPhasesSys, false);

    // Mark stable solution phases
    for (int i = 0; i < thermo.nElements; ++i) {
        int idx = thermo.iAssemblage(i);
        if (idx < 0) {
            // Negative index indicates solution phase
            int phaseIdx = -idx - 1;
            if (phaseIdx >= 0 && phaseIdx < static_cast<int>(gem.lSolnPhases.size())) {
                gem.lSolnPhases[phaseIdx] = true;
            }
        }
    }
}

void GEMSolver::checkSysOnlyPureConPhases(ThermoContext& ctx) {
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;
    auto& io = *ctx.io;

    // If only pure condensed phases, check mass balance
    if (thermo.nSolnPhases == 0 && thermo.nConPhases > 0) {
        // Calculate mass balance residuals
        double maxResidual = 0.0;

        for (int j = 0; j < thermo.nElements; ++j) {
            double sum = 0.0;
            for (int i = 0; i < thermo.nConPhases; ++i) {
                int speciesIdx = thermo.iAssemblage(i);
                if (speciesIdx > 0 && speciesIdx <= thermo.nSpecies) {
                    sum += thermo.dMolesPhase(i) * thermo.dStoichSpecies(speciesIdx - 1, j);
                }
            }

            double residual = std::abs(sum - thermo.dMolesElement(j));
            if (thermo.dMolesElement(j) > 1.0) {
                residual /= thermo.dMolesElement(j);
            }
            maxResidual = std::max(maxResidual, residual);
        }

        // Check if converged
        if (maxResidual < thermo.tolerances[kTolMassBalance]) {
            gem.lConverged = true;
        }
    }

    // If no solution phases yet and we have solution phases in the system,
    // add them if needed (e.g., gas phase for gas-dominated systems)
    if (thermo.nSolnPhases == 0 && thermo.nSolnPhasesSys > 0) {
        double R = Constants::kIdealGasConstant;
        double T = io.dTemperature;

        // Check if element potentials are valid (set by leveling)
        bool haveValidPotentials = (thermo.nConPhases > 0);
        if (haveValidPotentials) {
            for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
                if (thermo.dMolesElement(j) > 0 &&
                    (std::isnan(thermo.dElementPotential(j)) ||
                     std::isinf(thermo.dElementPotential(j)) ||
                     std::abs(thermo.dElementPotential(j)) > 1e100)) {
                    haveValidPotentials = false;
                    break;
                }
            }
        }

        // Try each solution phase
        for (int iPhase = 0; iPhase < thermo.nSolnPhasesSys; ++iPhase) {
            int iFirst = (iPhase > 0) ? thermo.nSpeciesPhase(iPhase) : 0;
            int iLast = thermo.nSpeciesPhase(iPhase + 1);
            int nPhaseSpecies = iLast - iFirst;
            if (nPhaseSpecies <= 0) continue;

            // Build list of feasible species (only contain active elements)
            std::vector<int> feasibleSpecies;
            for (int i = iFirst; i < iLast; ++i) {
                bool feasible = true;
                for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
                    if (thermo.dStoichSpecies(i, j) > 0.0 && thermo.dMolesElement(j) <= 0.0) {
                        feasible = false;
                        break;
                    }
                }
                if (feasible) {
                    feasibleSpecies.push_back(i);
                }
            }

            if (feasibleSpecies.empty()) continue;

            bool shouldAdd = false;
            std::vector<double> moleFractions(nPhaseSpecies, 0.0);

            if (haveValidPotentials) {
                // Use element potentials to estimate mole fractions
                // For ideal mixing: x_i ∝ exp(mu* - mu_std)
                double sumExp = 0.0;
                std::vector<double> expVal(nPhaseSpecies, 0.0);

                for (int i : feasibleSpecies) {
                    double muStar = 0.0;
                    for (int j = 0; j < thermo.nElements; ++j) {
                        muStar += thermo.dElementPotential(j) * thermo.dStoichSpecies(i, j);
                    }
                    muStar /= static_cast<double>(thermo.iParticlesPerMole(i));

                    double muStd = thermo.dStdGibbsEnergy(i) / (R * T);
                    double arg = muStar - muStd;
                    arg = std::max(-100.0, std::min(100.0, arg));
                    expVal[i - iFirst] = std::exp(arg);
                    sumExp += expVal[i - iFirst];
                }

                if (sumExp > 1e-300) {
                    // Compute driving force
                    double drivingForce = 0.0;
                    for (int i : feasibleSpecies) {
                        double x = expVal[i - iFirst] / sumExp;
                        moleFractions[i - iFirst] = x;
                        if (x > 1e-300) {
                            double muStar = 0.0;
                            for (int j = 0; j < thermo.nElements; ++j) {
                                muStar += thermo.dElementPotential(j) * thermo.dStoichSpecies(i, j);
                            }
                            muStar /= static_cast<double>(thermo.iParticlesPerMole(i));

                            double muStd = thermo.dStdGibbsEnergy(i) / (R * T);
                            double mu = muStd + std::log(x);
                            drivingForce += x * (mu - muStar);
                        }
                    }

                    // If driving force is favorable, add the phase
                    shouldAdd = (-drivingForce > thermo.tolerances[kTolDrivingForce]);
                }
            } else {
                // No valid element potentials - use stoichiometry-based initialization
                // Add the phase if it can satisfy mass balance for at least one element
                // Use Gibbs energy to weight species: x_i ∝ exp(-G_std_i / RT)
                double sumExp = 0.0;
                std::vector<double> expVal(nPhaseSpecies, 0.0);

                for (int i : feasibleSpecies) {
                    double muStd = thermo.dStdGibbsEnergy(i) / (R * T);
                    // Weight by negative Gibbs energy (lower is more favorable)
                    double arg = -muStd;
                    arg = std::max(-100.0, std::min(100.0, arg));
                    expVal[i - iFirst] = std::exp(arg);
                    sumExp += expVal[i - iFirst];
                }

                if (sumExp > 1e-300) {
                    for (int i : feasibleSpecies) {
                        moleFractions[i - iFirst] = expVal[i - iFirst] / sumExp;
                    }
                    // Add phase if it has feasible species
                    shouldAdd = true;
                }
            }

            if (shouldAdd) {
                // Add solution phase to assemblage
                PhaseAssemblage::addSolnPhase(ctx, iPhase);
                gem.iterLast = 0;

                // Initialize mole fractions
                for (int i = iFirst; i < iLast; ++i) {
                    thermo.dMolFraction(i) = moleFractions[i - iFirst];
                }

                // Estimate initial phase moles from mass balance
                double sumMoles = 0.0;
                int nActiveElements = 0;
                for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
                    if (thermo.dMolesElement(j) <= 0) continue;
                    nActiveElements++;

                    double avgStoich = 0.0;
                    for (int i : feasibleSpecies) {
                        avgStoich += thermo.dMolFraction(i) * thermo.dStoichSpecies(i, j);
                    }
                    if (avgStoich > 1e-20) {
                        sumMoles += thermo.dMolesElement(j) / avgStoich;
                    }
                }
                if (nActiveElements > 0) {
                    int assembIdx = thermo.nElements - thermo.nSolnPhases;
                    if (assembIdx >= 0 && assembIdx < thermo.dMolesPhase.size()) {
                        thermo.dMolesPhase(assembIdx) = sumMoles / nActiveElements;
                    }
                }

                break;  // Only add one phase at a time
            }
        }
    }
}

} // namespace Thermochimica
