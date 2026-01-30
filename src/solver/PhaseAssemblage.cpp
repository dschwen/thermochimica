#include "thermochimica/solver/GEMSolver.hpp"
#include "thermochimica/util/Constants.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <limits>

namespace Thermochimica {

/// @brief Check if a species is feasible (contains only active elements)
static bool speciesIsFeasibleForDF(ThermoContext& ctx, int speciesIdx) {
    auto& thermo = *ctx.thermo;
    int nEl = thermo.nElements - thermo.nChargedConstraints;

    for (int j = 0; j < nEl; ++j) {
        if (thermo.dStoichSpecies(speciesIdx, j) > 0.0 && thermo.dMolesElement(j) <= 0.0) {
            return false;  // Species contains an element not in the input
        }
    }
    return true;
}

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

        // Build list of feasible species for this phase
        std::vector<int> feasibleSpecies;
        for (int i = iFirst; i < iLast; ++i) {
            if (speciesIsFeasibleForDF(ctx, i)) {
                feasibleSpecies.push_back(i);
            }
        }

        if (feasibleSpecies.empty()) {
            gem.dDrivingForceSoln(iPhase) = 0.0;  // No feasible species
            continue;
        }

        // For each feasible species, compute mu* from element potentials
        // and estimate optimal mole fractions
        double sumExp = 0.0;
        std::vector<double> expVal(iLast - iFirst, 0.0);

        for (int i : feasibleSpecies) {
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
            for (int i : feasibleSpecies) {
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

    // Check phase stability - allow Newton solver to converge before adding phases
    // This prevents premature phase additions with incorrect element potentials
    // which causes oscillation between phases
    int checkInterval = (gem.iterGlobal < 100) ? 50 : gem.iterStep * 5;
    if (gem.iterGlobal - gem.iterLast < checkInterval) {
        return;
    }

    // Compute driving forces for all unstable solution phases
    computeSolutionPhaseDrivingForces(ctx);

    // Find phase with highest driving force (most favorable to add)
    int bestPhase = -1;
    double bestDF = thermo.tolerances[kTolDrivingForce];
    for (int iPhase = 0; iPhase < thermo.nSolnPhasesSys; ++iPhase) {
        if (gem.lSolnPhases[iPhase]) {
            continue;  // Already stable
        }

        double drivingForce = gem.dDrivingForceSoln(iPhase);
        if (drivingForce > bestDF) {
            bestDF = drivingForce;
            bestPhase = iPhase;
        }
    }

    // Try to add the best phase
    if (bestPhase >= 0) {
        // Try to add the solution phase
            if (addSolnPhase(ctx, bestPhase)) {
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
                if (worstPhaseIdx >= 0 && bestDF > -worstDF) {
                    int removedSpecies = thermo.iAssemblage(worstPhaseIdx) - 1;
                    removePureConPhase(ctx, removedSpecies);

                    if (addSolnPhase(ctx, bestPhase)) {
                        gem.iterLast = gem.iterGlobal;
                        gem.iterLastSoln = gem.iterGlobal;
                        return;
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

/// @brief Check if a species is feasible (contains only active elements)
static bool speciesIsFeasibleForPhase(ThermoContext& ctx, int speciesIdx) {
    auto& thermo = *ctx.thermo;
    int nEl = thermo.nElements - thermo.nChargedConstraints;

    for (int j = 0; j < nEl; ++j) {
        if (thermo.dStoichSpecies(speciesIdx, j) > 0.0 && thermo.dMolesElement(j) <= 0.0) {
            return false;  // Species contains an element not in the input
        }
    }
    return true;
}

bool PhaseAssemblage::addSolnPhase(ThermoContext& ctx, int phaseIndex) {
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;

    // Count active elements (with non-zero input)
    int nActiveElements = 0;
    for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
        if (thermo.dMolesElement(j) > 0.0) {
            ++nActiveElements;
        }
    }

    // Check phase rule: max phases = active components
    if (thermo.nConPhases + thermo.nSolnPhases >= nActiveElements) {
        return false;  // Would violate phase rule
    }

    // Add phase to assemblage
    thermo.nSolnPhases++;
    int newIdx = thermo.nElements - thermo.nSolnPhases;
    thermo.iAssemblage(newIdx) = -(phaseIndex + 1);

    // Get species range for this phase
    int iFirst = (phaseIndex > 0) ? thermo.nSpeciesPhase(phaseIndex) : 0;
    int iLast = thermo.nSpeciesPhase(phaseIndex + 1);

    // Build list of feasible species (only contain active elements)
    std::vector<int> feasibleSpecies;
    for (int i = iFirst; i < iLast; ++i) {
        if (speciesIsFeasibleForPhase(ctx, i)) {
            feasibleSpecies.push_back(i);
        }
    }

    if (feasibleSpecies.empty()) {
        // No feasible species - remove the phase we just added
        thermo.nSolnPhases--;
        thermo.iAssemblage(newIdx) = 0;
        return false;
    }

    // Estimate initial moles based on mass balance requirements
    // First, compute what the existing phases (pure + solution) provide
    Eigen::VectorXd accountedMass = Eigen::VectorXd::Zero(thermo.nElements);

    // Pure condensed phases
    for (int i = 0; i < thermo.nConPhases; ++i) {
        int speciesIdx = thermo.iAssemblage(i) - 1;
        if (speciesIdx >= 0 && speciesIdx < thermo.nSpecies) {
            for (int j = 0; j < thermo.nElements; ++j) {
                accountedMass(j) += thermo.dMolesPhase(i) * thermo.dStoichSpecies(speciesIdx, j);
            }
        }
    }

    // Already-stable solution phases (nSolnPhases-1 because we just incremented for the new phase)
    for (int iPhase = 0; iPhase < thermo.nSolnPhases - 1; ++iPhase) {
        int assembIdx = thermo.nElements - (thermo.nSolnPhases - 1) + iPhase;
        int solnIdx = -thermo.iAssemblage(assembIdx) - 1;
        if (solnIdx < 0 || solnIdx >= thermo.nSolnPhasesSys) continue;

        int iFirst = (solnIdx > 0) ? thermo.nSpeciesPhase(solnIdx) : 0;
        int iLast = thermo.nSpeciesPhase(solnIdx + 1);
        double phaseMoles = thermo.dMolesPhase(assembIdx);

        for (int i = iFirst; i < iLast; ++i) {
            if (thermo.dMolFraction(i) > 0) {
                for (int j = 0; j < thermo.nElements; ++j) {
                    accountedMass(j) += phaseMoles * thermo.dMolFraction(i) * thermo.dStoichSpecies(i, j);
                }
            }
        }
    }

    // Compute deficit for each element
    Eigen::VectorXd deficit = Eigen::VectorXd::Zero(thermo.nElements);
    double totalDeficit = 0.0;
    for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
        if (thermo.dMolesElement(j) > 0.0) {
            deficit(j) = std::max(0.0, thermo.dMolesElement(j) - accountedMass(j));
            totalDeficit += deficit(j);
        }
    }

    // If there's no mass deficit, the new phase isn't needed
    if (totalDeficit < thermo.tolerances[kTolPhaseMoles]) {
        // Phase not needed - remove it
        thermo.nSolnPhases--;
        thermo.iAssemblage(newIdx) = 0;
        gem.lSolnPhases[phaseIndex] = false;
        // Also unmark miscibility partners
        for (int p = 0; p < thermo.nSolnPhasesSys; ++p) {
            if (thermo.cSolnPhaseName[p] == thermo.cSolnPhaseName[phaseIndex]) {
                gem.lSolnPhases[p] = false;
            }
        }
        return false;
    }

    // Estimate moles needed based on expected dominant species
    double initialMoles = 0.1;

    // Find the feasible species with lowest chemical potential (most stable)
    int bestSpecies = feasibleSpecies[0];
    double bestMu = std::numeric_limits<double>::max();
    for (int i : feasibleSpecies) {
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
    // Use overall composition as starting point (not element potentials,
    // which can give degenerate solutions like pure end-members)

    // First, compute total active element moles
    double totalActiveMoles = 0.0;
    for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
        if (thermo.dMolesElement(j) > 0.0) {
            totalActiveMoles += thermo.dMolesElement(j);
        }
    }

    // Initialize mole fractions based on system composition
    // For each feasible species, estimate its contribution based on
    // how well it matches the overall composition
    std::vector<double> moleFracs(iLast - iFirst, 1e-100);
    double sumFrac = 0.0;

    for (int i : feasibleSpecies) {
        // Sum the stoichiometry for active elements
        double stoichSum = 0.0;
        for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
            if (thermo.dMolesElement(j) > 0.0) {
                stoichSum += thermo.dStoichSpecies(i, j);
            }
        }

        if (stoichSum > 0) {
            // Weight by how much of the system this species could represent
            double weight = 0.0;
            for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
                if (thermo.dMolesElement(j) > 0.0 && thermo.dStoichSpecies(i, j) > 0) {
                    // Mole fraction based on element j
                    weight += (thermo.dMolesElement(j) / totalActiveMoles) *
                              (thermo.dStoichSpecies(i, j) / stoichSum);
                }
            }
            moleFracs[i - iFirst] = std::max(weight, 1e-10);
            sumFrac += moleFracs[i - iFirst];
        }
    }

    // Normalize mole fractions
    if (sumFrac > 1e-300) {
        for (int i = iFirst; i < iLast; ++i) {
            thermo.dMolFraction(i) = moleFracs[i - iFirst] / sumFrac;
            thermo.dMolesSpecies(i) = initialMoles * thermo.dMolFraction(i);
        }
    } else {
        // Fallback: equal mole fractions among feasible species
        int nFeasible = static_cast<int>(feasibleSpecies.size());
        for (int i = iFirst; i < iLast; ++i) {
            thermo.dMolFraction(i) = 1e-100;
            thermo.dMolesSpecies(i) = initialMoles * 1e-100;
        }
        if (nFeasible > 0) {
            double x = 1.0 / nFeasible;
            for (int i : feasibleSpecies) {
                thermo.dMolFraction(i) = x;
                thermo.dMolesSpecies(i) = initialMoles * x;
            }
        }
    }

    // Mark phase as stable
    gem.lSolnPhases[phaseIndex] = true;

    // Also mark the paired miscibility phase (if any) as stable
    // This prevents the solver from trying to add both miscibility partners
    // which would cause one to have zero moles and get removed
    for (int p = 0; p < thermo.nSolnPhasesSys; ++p) {
        if (p != phaseIndex &&
            thermo.cSolnPhaseName[p] == thermo.cSolnPhaseName[phaseIndex]) {
            gem.lSolnPhases[p] = true;
        }
    }

    // Now solve for phase moles that satisfy mass balance
    // Build stoichiometry matrix for current assemblage (pure + solution phases)
    int nPhasesTotal = thermo.nConPhases + 1;  // +1 for the just-added solution phase

    // Count active elements and build list
    std::vector<int> activeElementList;
    for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
        if (thermo.dMolesElement(j) > 0.0) {
            activeElementList.push_back(j);
        }
    }
    int nActiveEl = static_cast<int>(activeElementList.size());

    if (nPhasesTotal == nActiveEl) {
        // Square system - can solve directly
        Eigen::MatrixXd stoichMat = Eigen::MatrixXd::Zero(nActiveEl, nPhasesTotal);
        Eigen::VectorXd targetMoles = Eigen::VectorXd::Zero(nActiveEl);

        // Pure condensed phases
        for (int p = 0; p < thermo.nConPhases; ++p) {
            int speciesIdx = thermo.iAssemblage(p) - 1;
            if (speciesIdx >= 0 && speciesIdx < thermo.nSpecies) {
                for (int jRed = 0; jRed < nActiveEl; ++jRed) {
                    int jFull = activeElementList[jRed];
                    stoichMat(jRed, p) = thermo.dStoichSpecies(speciesIdx, jFull);
                }
            }
        }

        // Solution phase (use average stoichiometry weighted by mole fractions)
        for (int jRed = 0; jRed < nActiveEl; ++jRed) {
            int jFull = activeElementList[jRed];
            double avgStoich = 0.0;
            for (int i = iFirst; i < iLast; ++i) {
                if (thermo.dMolFraction(i) > 0) {
                    avgStoich += thermo.dMolFraction(i) * thermo.dStoichSpecies(i, jFull);
                }
            }
            stoichMat(jRed, thermo.nConPhases) = avgStoich;
        }

        for (int jRed = 0; jRed < nActiveEl; ++jRed) {
            int jFull = activeElementList[jRed];
            targetMoles(jRed) = thermo.dMolesElement(jFull);
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

        // Check if the new phase would have significant moles
        // If the stoichiometry solve gives zero moles, the phase isn't needed
        double newPhaseMoles = phaseMoles(thermo.nConPhases);
        if (newPhaseMoles < thermo.tolerances[kTolPhaseMoles] * 10.0) {
            // Phase not needed - remove it
            thermo.nSolnPhases--;
            thermo.iAssemblage(newIdx) = 0;
            gem.lSolnPhases[phaseIndex] = false;
            // Also unmark miscibility partners
            for (int p = 0; p < thermo.nSolnPhasesSys; ++p) {
                if (thermo.cSolnPhaseName[p] == thermo.cSolnPhaseName[phaseIndex]) {
                    gem.lSolnPhases[p] = false;
                }
            }
            return false;
        }

        thermo.dMolesPhase(newIdx) = std::max(1e-10, newPhaseMoles);

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

    // Zero out species moles and mole fractions for this phase
    // (important for mass balance and phase fraction calculations!)
    int iFirst = (phaseIndex > 0) ? thermo.nSpeciesPhase(phaseIndex) : 0;
    int iLast = thermo.nSpeciesPhase(phaseIndex + 1);
    for (int i = iFirst; i < iLast; ++i) {
        thermo.dMolesSpecies(i) = 0.0;
        thermo.dMolFraction(i) = 0.0;
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

    // Count active elements (with non-zero input)
    int nActiveElements = 0;
    for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
        if (thermo.dMolesElement(j) > 0.0) {
            ++nActiveElements;
        }
    }

    // Check phase rule: max phases = active components
    if (thermo.nConPhases + thermo.nSolnPhases >= nActiveElements) {
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
