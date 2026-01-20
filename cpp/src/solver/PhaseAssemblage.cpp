#include "thermochimica/solver/GEMSolver.hpp"
#include "thermochimica/util/Constants.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <limits>

namespace Thermochimica {

/// @brief Compute driving forces for all unstable solution phases
/// @details For each solution phase not currently in the assemblage,
/// estimates optimal mole fractions and computes the driving force
/// to determine if the phase should be added.
static void computeSolutionPhaseDrivingForces(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;
    auto& gem = *ctx.gem;

    // Compute driving force for each solution phase not currently stable
    for (int iPhase = 0; iPhase < thermo.nSolnPhasesSys; ++iPhase) {
        if (gem.lSolnPhases[iPhase]) {
            gem.dDrivingForceSoln(iPhase) = 0.0;  // Already stable
            continue;
        }

        // Get species range for this phase
        int iFirst = (iPhase > 0) ? thermo.nSpeciesPhase(iPhase) : 0;
        int iLast = thermo.nSpeciesPhase(iPhase + 1);

        if (iFirst >= iLast) {
            gem.dDrivingForceSoln(iPhase) = 0.0;
            continue;
        }

        // For each species, compute mu* from element potentials
        // and estimate optimal mole fractions
        double sumExp = 0.0;
        std::vector<double> expVal(iLast - iFirst, 0.0);

        for (int i = iFirst; i < iLast; ++i) {
            // Compute chemical potential from element potentials
            // mu*_i = sum_j (lambda_j * a_ij) / n_i
            double muStar = 0.0;
            for (int j = 0; j < thermo.nElements; ++j) {
                muStar += thermo.dElementPotential(j) * thermo.dStoichSpecies(i, j);
            }
            muStar /= static_cast<double>(thermo.iParticlesPerMole(i));

            // Standard Gibbs energy (already normalized by RT)
            double muStd = thermo.dStdGibbsEnergy(i) /
                          (Constants::kIdealGasConstant * ctx.io->dTemperature);

            // Estimate mole fraction: x_i ∝ exp(mu* - mu_std)
            double arg = muStar - muStd;
            arg = std::max(-100.0, std::min(100.0, arg));  // Prevent overflow
            expVal[i - iFirst] = std::exp(arg);
            sumExp += expVal[i - iFirst];
        }

        // Normalize and compute driving force
        // DF = sum_i x_i * (mu_std + ln(x_i) - mu*)
        double drivingForce = 0.0;

        if (sumExp > 1e-300) {
            for (int i = iFirst; i < iLast; ++i) {
                double x = expVal[i - iFirst] / sumExp;
                if (x > 1e-300) {
                    // Compute chemical potential from element potentials
                    double muStar = 0.0;
                    for (int j = 0; j < thermo.nElements; ++j) {
                        muStar += thermo.dElementPotential(j) * thermo.dStoichSpecies(i, j);
                    }
                    muStar /= static_cast<double>(thermo.iParticlesPerMole(i));

                    // Standard Gibbs energy (normalized by RT)
                    double muStd = thermo.dStdGibbsEnergy(i) /
                                  (Constants::kIdealGasConstant * ctx.io->dTemperature);

                    // Chemical potential = mu_std + ln(x) for ideal mixing
                    double mu = muStd + std::log(x);

                    // Driving force contribution: x_i * (mu_i - mu*_i)
                    // Negative driving force means phase is favorable
                    drivingForce += x * (mu - muStar);
                }
            }
        }

        // Store driving force (negative = favorable)
        // The check function looks for drivingForce > tolerance, so we negate
        gem.dDrivingForceSoln(iPhase) = -drivingForce;

    }
}

void PhaseAssemblage::check(ThermoContext& ctx) {
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;

    // Check if system should be reverted
    if (gem.lRevertSystem) {
        revert(ctx);
        gem.lRevertSystem = false;
        return;
    }

    // Check for phases that should be removed (negative or very small moles)
    for (int i = 0; i < thermo.nSolnPhases; ++i) {
        int phaseIdx = thermo.nElements - thermo.nSolnPhases + i;
        if (thermo.dMolesPhase(phaseIdx) < thermo.tolerances[kTolPhaseMoles]) {
            int solnIdx = -thermo.iAssemblage(phaseIdx) - 1;
            if (solnIdx >= 0) {
                removeSolnPhase(ctx, solnIdx);
                gem.iterLast = gem.iterGlobal;
                return;
            }
        }
    }

    for (int i = 0; i < thermo.nConPhases; ++i) {
        if (thermo.dMolesPhase(i) < 0 ||
            thermo.dMolesPhase(i) < thermo.tolerances[kTolPhaseMoles]) {
            int speciesIdx = thermo.iAssemblage(i);
            if (speciesIdx > 0) {
                removePureConPhase(ctx, speciesIdx - 1);
                gem.iterLast = gem.iterGlobal;
                return;
            }
        }
    }

    // Check for phases that should be added (positive driving force)
    // Set default step size if not set
    if (gem.iterStep == 0) {
        gem.iterStep = 10;
    }

    // Only check periodically to avoid oscillation
    if (gem.iterGlobal - gem.iterLast < gem.iterStep * 2) {
        return;
    }

    // Compute driving forces for all unstable solution phases
    computeSolutionPhaseDrivingForces(ctx);

    // Check solution phases
    for (int iPhase = 0; iPhase < thermo.nSolnPhasesSys; ++iPhase) {
        if (gem.lSolnPhases[iPhase]) {
            continue;  // Already stable
        }

        double drivingForce = gem.dDrivingForceSoln(iPhase);
        if (drivingForce > thermo.tolerances[kTolDrivingForce]) {
            // Try to add the solution phase
            if (addSolnPhase(ctx, iPhase)) {
                gem.iterLast = gem.iterGlobal;
                gem.iterLastSoln = gem.iterGlobal;
                return;
            }

            // If we couldn't add due to phase rule, try swapping with worst pure phase
            if (thermo.nConPhases > 0) {
                // Find the pure condensed phase with worst (most positive) driving force
                int worstPhaseIdx = -1;
                double worstDF = -1e30;

                for (int i = 0; i < thermo.nConPhases; ++i) {
                    int speciesIdx = thermo.iAssemblage(i) - 1;  // Convert to 0-based
                    if (speciesIdx < 0 || speciesIdx >= thermo.nSpecies) continue;

                    // Compute driving force for this pure phase
                    double muCalc = 0.0;
                    for (int j = 0; j < thermo.nElements; ++j) {
                        muCalc += thermo.dElementPotential(j) * thermo.dStoichSpecies(speciesIdx, j);
                    }
                    muCalc /= thermo.iParticlesPerMole(speciesIdx);

                    double phaseDF = thermo.dStdGibbsEnergy(speciesIdx) - muCalc;

                    // Positive driving force means phase is unstable
                    if (phaseDF > worstDF) {
                        worstDF = phaseDF;
                        worstPhaseIdx = i;
                    }
                }

                // If solution phase has better driving force than worst pure phase, swap
                if (worstPhaseIdx >= 0 && drivingForce > -worstDF) {
                    int removedSpecies = thermo.iAssemblage(worstPhaseIdx) - 1;
                    removePureConPhase(ctx, removedSpecies);

                    if (addSolnPhase(ctx, iPhase)) {
                        gem.iterLast = gem.iterGlobal;
                        gem.iterLastSoln = gem.iterGlobal;
                        return;
                    }
                }
            }
        }
    }

    // Check pure condensed phases
    for (int i = thermo.nSpeciesPhase(thermo.nSolnPhasesSys); i < thermo.nSpecies; ++i) {
        // Calculate driving force for pure phase
        double muCalc = 0.0;
        for (int j = 0; j < thermo.nElements; ++j) {
            muCalc += thermo.dElementPotential(j) * thermo.dStoichSpecies(i, j);
        }
        muCalc /= thermo.iParticlesPerMole(i);

        double drivingForce = thermo.dStdGibbsEnergy(i) - muCalc;

        if (drivingForce < -thermo.tolerances[kTolDrivingForce]) {
            // Phase is more stable than current assemblage
            if (addPureConPhase(ctx, i)) {
                gem.iterLast = gem.iterGlobal;
                gem.iterLastCon = gem.iterGlobal;
                return;
            }
        }
    }
}

bool PhaseAssemblage::addSolnPhase(ThermoContext& ctx, int phaseIndex) {
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;

    // Check phase rule
    if (thermo.nConPhases + thermo.nSolnPhases >= thermo.nElements - thermo.nChargedConstraints) {
        return false;  // Would violate phase rule
    }

    // Add phase to assemblage
    thermo.nSolnPhases++;
    int newIdx = thermo.nElements - thermo.nSolnPhases;
    thermo.iAssemblage(newIdx) = -(phaseIndex + 1);

    // Estimate initial moles based on mass balance requirements
    // First, compute what the pure phases provide
    Eigen::VectorXd accountedMass = Eigen::VectorXd::Zero(thermo.nElements);
    for (int i = 0; i < thermo.nConPhases; ++i) {
        int speciesIdx = thermo.iAssemblage(i) - 1;
        if (speciesIdx >= 0 && speciesIdx < thermo.nSpecies) {
            for (int j = 0; j < thermo.nElements; ++j) {
                accountedMass(j) += thermo.dMolesPhase(i) * thermo.dStoichSpecies(speciesIdx, j);
            }
        }
    }

    // Compute deficit for each element
    Eigen::VectorXd deficit = Eigen::VectorXd::Zero(thermo.nElements);
    for (int j = 0; j < thermo.nElements; ++j) {
        deficit(j) = std::max(0.0, thermo.dMolesElement(j) - accountedMass(j));
    }

    // Estimate moles needed based on expected dominant species (e.g., CO2)
    // Use a simple heuristic: divide by average stoichiometry of deficient elements
    double initialMoles = 0.1;
    int iFirst = (phaseIndex > 0) ? thermo.nSpeciesPhase(phaseIndex) : 0;
    int iLast = thermo.nSpeciesPhase(phaseIndex + 1);

    // Find the species with lowest chemical potential (most stable)
    int bestSpecies = iFirst;
    double bestMu = std::numeric_limits<double>::max();
    for (int i = iFirst; i < iLast; ++i) {
        double muStd = thermo.dStdGibbsEnergy(i) /
                      (Constants::kIdealGasConstant * ctx.io->dTemperature);
        if (muStd < bestMu) {
            bestMu = muStd;
            bestSpecies = i;
        }
    }

    // Estimate moles based on most stable species' stoichiometry
    for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
        if (deficit(j) > 0 && thermo.dStoichSpecies(bestSpecies, j) > 0) {
            double molesNeeded = deficit(j) / thermo.dStoichSpecies(bestSpecies, j);
            initialMoles = std::max(initialMoles, molesNeeded);
        }
    }

    thermo.dMolesPhase(newIdx) = initialMoles;

    // Initialize mole fractions for species in this phase
    // Use optimal mole fractions from element potentials
    // (iFirst and iLast already computed above)

    double sumExp = 0.0;
    std::vector<double> expVal(iLast - iFirst, 0.0);

    for (int i = iFirst; i < iLast; ++i) {
        // Compute mu* from element potentials
        double muStar = 0.0;
        for (int j = 0; j < thermo.nElements; ++j) {
            muStar += thermo.dElementPotential(j) * thermo.dStoichSpecies(i, j);
        }
        muStar /= static_cast<double>(thermo.iParticlesPerMole(i));

        // Standard Gibbs energy (normalized by RT)
        double muStd = thermo.dStdGibbsEnergy(i) /
                      (Constants::kIdealGasConstant * ctx.io->dTemperature);

        // Optimal mole fraction estimate
        double arg = muStar - muStd;
        arg = std::max(-100.0, std::min(100.0, arg));
        expVal[i - iFirst] = std::exp(arg);
        sumExp += expVal[i - iFirst];
    }

    // Normalize mole fractions and set species moles
    if (sumExp > 1e-300) {
        for (int i = iFirst; i < iLast; ++i) {
            double x = expVal[i - iFirst] / sumExp;
            x = std::max(x, 1e-100);  // Minimum mole fraction
            thermo.dMolFraction(i) = x;
            thermo.dMolesSpecies(i) = initialMoles * x;
        }
    } else {
        // Equal mole fractions as fallback
        int nSpecies = iLast - iFirst;
        for (int i = iFirst; i < iLast; ++i) {
            thermo.dMolFraction(i) = 1.0 / nSpecies;
            thermo.dMolesSpecies(i) = initialMoles / nSpecies;
        }
    }

    // Mark phase as stable
    gem.lSolnPhases[phaseIndex] = true;

    // Now solve for phase moles that satisfy mass balance
    // Build stoichiometry matrix for current assemblage (pure + solution phases)
    int nPhases = thermo.nConPhases + 1;  // +1 for the just-added solution phase
    int nEl = thermo.nElements - thermo.nChargedConstraints;

    if (nPhases == nEl) {
        // Square system - can solve directly
        Eigen::MatrixXd stoichMat = Eigen::MatrixXd::Zero(nEl, nPhases);
        Eigen::VectorXd targetMoles = Eigen::VectorXd::Zero(nEl);

        // Pure condensed phases
        for (int p = 0; p < thermo.nConPhases; ++p) {
            int speciesIdx = thermo.iAssemblage(p) - 1;
            if (speciesIdx >= 0 && speciesIdx < thermo.nSpecies) {
                for (int j = 0; j < nEl; ++j) {
                    stoichMat(j, p) = thermo.dStoichSpecies(speciesIdx, j);
                }
            }
        }

        // Solution phase (use average stoichiometry weighted by mole fractions)
        for (int j = 0; j < nEl; ++j) {
            double avgStoich = 0.0;
            for (int i = iFirst; i < iLast; ++i) {
                if (thermo.dMolFraction(i) > 0) {
                    avgStoich += thermo.dMolFraction(i) * thermo.dStoichSpecies(i, j);
                }
            }
            stoichMat(j, thermo.nConPhases) = avgStoich;
        }

        for (int j = 0; j < nEl; ++j) {
            targetMoles(j) = thermo.dMolesElement(j);
        }

        // Solve stoichMat @ phaseMoles = targetMoles
        Eigen::VectorXd phaseMoles = stoichMat.colPivHouseholderQr().solve(targetMoles);

        // Apply the solution (clamp to non-negative)
        for (int p = 0; p < thermo.nConPhases; ++p) {
            thermo.dMolesPhase(p) = std::max(0.0, phaseMoles(p));
            int speciesIdx = thermo.iAssemblage(p) - 1;
            if (speciesIdx >= 0) {
                thermo.dMolesSpecies(speciesIdx) = thermo.dMolesPhase(p);
            }
        }
        thermo.dMolesPhase(newIdx) = std::max(1e-10, phaseMoles(thermo.nConPhases));

        // Update species moles for solution phase
        for (int i = iFirst; i < iLast; ++i) {
            if (thermo.dMolFraction(i) > 0) {
                thermo.dMolesSpecies(i) = thermo.dMolesPhase(newIdx) * thermo.dMolFraction(i);
            }
        }
    }

    // Store in history
    for (int i = 0; i < thermo.nElements; ++i) {
        gem.iterHistory(i, gem.iterGlobal) = thermo.iAssemblage(i);
    }

    return true;
}

bool PhaseAssemblage::removeSolnPhase(ThermoContext& ctx, int phaseIndex) {
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;

    if (thermo.nSolnPhases <= 0) {
        return false;
    }

    // Find phase in assemblage
    int removeIdx = -1;
    for (int i = thermo.nElements - thermo.nSolnPhases; i < thermo.nElements; ++i) {
        if (thermo.iAssemblage(i) == -(phaseIndex + 1)) {
            removeIdx = i;
            break;
        }
    }

    if (removeIdx < 0) {
        return false;
    }

    // Shift remaining phases
    for (int i = removeIdx; i < thermo.nElements - 1; ++i) {
        thermo.iAssemblage(i) = thermo.iAssemblage(i + 1);
        thermo.dMolesPhase(i) = thermo.dMolesPhase(i + 1);
    }

    thermo.iAssemblage(thermo.nElements - 1) = 0;
    thermo.dMolesPhase(thermo.nElements - 1) = 0.0;
    thermo.nSolnPhases--;

    // Mark phase as unstable
    gem.lSolnPhases[phaseIndex] = false;

    return true;
}

bool PhaseAssemblage::addPureConPhase(ThermoContext& ctx, int speciesIndex) {
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;

    // Check phase rule
    if (thermo.nConPhases + thermo.nSolnPhases >= thermo.nElements - thermo.nChargedConstraints) {
        return false;
    }

    // Check if species already in assemblage
    for (int i = 0; i < thermo.nConPhases; ++i) {
        if (thermo.iAssemblage(i) == speciesIndex + 1) {
            return false;  // Already present
        }
    }

    // Add to assemblage
    thermo.iAssemblage(thermo.nConPhases) = speciesIndex + 1;
    thermo.dMolesPhase(thermo.nConPhases) = thermo.tolerances[kTolPhaseMoles] * 10.0;
    thermo.nConPhases++;

    // Store in history
    for (int i = 0; i < thermo.nElements; ++i) {
        gem.iterHistory(i, gem.iterGlobal) = thermo.iAssemblage(i);
    }

    return true;
}

bool PhaseAssemblage::removePureConPhase(ThermoContext& ctx, int speciesIndex) {
    auto& thermo = *ctx.thermo;

    if (thermo.nConPhases <= 0) {
        return false;
    }

    // Find species in assemblage
    int removeIdx = -1;
    for (int i = 0; i < thermo.nConPhases; ++i) {
        if (thermo.iAssemblage(i) == speciesIndex + 1) {
            removeIdx = i;
            break;
        }
    }

    if (removeIdx < 0) {
        return false;
    }

    // Zero out the species moles (important for mass balance check!)
    thermo.dMolesSpecies(speciesIndex) = 0.0;
    thermo.dMolFraction(speciesIndex) = 0.0;

    // Shift remaining phases
    for (int i = removeIdx; i < thermo.nConPhases - 1; ++i) {
        thermo.iAssemblage(i) = thermo.iAssemblage(i + 1);
        thermo.dMolesPhase(i) = thermo.dMolesPhase(i + 1);
    }

    thermo.iAssemblage(thermo.nConPhases - 1) = 0;
    thermo.dMolesPhase(thermo.nConPhases - 1) = 0.0;
    thermo.nConPhases--;

    return true;
}

void PhaseAssemblage::swapPhases(ThermoContext& ctx) {
    // Placeholder for phase swapping logic
    // This handles cases where one phase should replace another
}

void PhaseAssemblage::revert(ThermoContext& ctx) {
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;

    // Revert to a previous successful assemblage
    int revertIter = std::max(1, gem.iterGlobal - 200);

    if (gem.iterGlobal > 1500) {
        revertIter = std::max(1, gem.iterGlobal - 1000);
    }

    // Restore assemblage from history
    for (int i = 0; i < thermo.nElements; ++i) {
        thermo.iAssemblage(i) = gem.iterHistory(i, revertIter);
    }

    // Recalculate phase counts
    thermo.nConPhases = 0;
    thermo.nSolnPhases = 0;

    for (int i = 0; i < thermo.nElements; ++i) {
        int idx = thermo.iAssemblage(i);
        if (idx > 0) {
            thermo.nConPhases++;
        } else if (idx < 0) {
            thermo.nSolnPhases++;
        }
    }

    // Update stable phase flags
    std::fill(gem.lSolnPhases.begin(), gem.lSolnPhases.end(), false);
    for (int i = 0; i < thermo.nElements; ++i) {
        int idx = thermo.iAssemblage(i);
        if (idx < 0) {
            int phaseIdx = -idx - 1;
            if (phaseIdx >= 0 && phaseIdx < static_cast<int>(gem.lSolnPhases.size())) {
                gem.lSolnPhases[phaseIdx] = true;
            }
        }
    }

    gem.iterRevert = gem.iterGlobal;
    gem.dGEMFunctionNorm *= 10.0;  // Relax convergence tolerance
}

} // namespace Thermochimica
