/// @file PhaseAssemblageManager.cpp
/// @brief Implementation of phase assemblage manager
/// @details Adapted from PhaseAssemblage.cpp static methods

#include "thermochimica/solver/PhaseAssemblageManager.hpp"
#include "thermochimica/context/ThermoState.hpp"
#include "thermochimica/context/GEMState.hpp"
#include "thermochimica/context/ThermoIO.hpp"
#include "thermochimica/util/Constants.hpp"
#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>

namespace Thermochimica {

PhaseAssemblageManager::PhaseAssemblageManager(ThermoState& state,
                                               GEMState& gemState,
                                               ThermoIO& io)
    : state_(state), gemState_(gemState), io_(io) {
}

bool PhaseAssemblageManager::speciesIsFeasible(int speciesIdx) const {
    int nEl = state_.nElements - state_.nChargedConstraints;

    for (int j = 0; j < nEl; ++j) {
        if (state_.dStoichSpecies(speciesIdx, j) > 0.0 &&
            state_.dMolesElement(j) <= 0.0) {
            return false;  // Species contains an element not in the input
        }
    }
    return true;
}

void PhaseAssemblageManager::computeSolutionPhaseDrivingForces() {
    // Compute driving force for each solution phase not currently stable
    for (int iPhase = 0; iPhase < state_.nSolnPhasesSys; ++iPhase) {
        if (gemState_.lSolnPhases[iPhase]) {
            gemState_.dDrivingForceSoln(iPhase) = 0.0;  // Already stable
            continue;
        }

        // Get species range for this phase
        int iFirst = (iPhase > 0) ? state_.nSpeciesPhase(iPhase) : 0;
        int iLast = state_.nSpeciesPhase(iPhase + 1);

        if (iFirst >= iLast) {
            gemState_.dDrivingForceSoln(iPhase) = 0.0;
            continue;
        }

        // Build list of feasible species for this phase
        std::vector<int> feasibleSpecies;
        for (int i = iFirst; i < iLast; ++i) {
            if (speciesIsFeasible(i)) {
                feasibleSpecies.push_back(i);
            }
        }

        if (feasibleSpecies.empty()) {
            gemState_.dDrivingForceSoln(iPhase) = 0.0;  // No feasible species
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
            for (int j = 0; j < state_.nElements; ++j) {
                muStar += state_.dElementPotential(j) * state_.dStoichSpecies(i, j);
            }
            muStar /= static_cast<double>(state_.iParticlesPerMole(i));

            // Standard Gibbs energy (already normalized by RT)
            double muStd = state_.dStdGibbsEnergy(i) /
                          (Constants::kIdealGasConstant * io_.dTemperature);

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
                    for (int j = 0; j < state_.nElements; ++j) {
                        muStar += state_.dElementPotential(j) * state_.dStoichSpecies(i, j);
                    }
                    muStar /= static_cast<double>(state_.iParticlesPerMole(i));

                    double muStd = state_.dStdGibbsEnergy(i) /
                                  (Constants::kIdealGasConstant * io_.dTemperature);

                    // Contribution to driving force (normalized by RT)
                    drivingForce += x * (muStd + std::log(x) - muStar);
                }
            }
        }

        gemState_.dDrivingForceSoln(iPhase) = -drivingForce;  // Negative for favorable
    }
}

void PhaseAssemblageManager::check() {
    // Check if system should be reverted
    if (gemState_.lRevertSystem) {
        revert();
        gemState_.lRevertSystem = false;
        return;
    }

    // Check for phases that should be removed (negative or very small moles)
    for (int i = 0; i < state_.nSolnPhases; ++i) {
        int phaseIdx = state_.nElements - state_.nSolnPhases + i;
        if (state_.dMolesPhase(phaseIdx) < state_.tolerances[kTolPhaseMoles]) {
            int solnIdx = -state_.iAssemblage(phaseIdx) - 1;
            if (solnIdx >= 0) {
                // Pass assemblage index i, not system phase index solnIdx
                removeSolnPhase(i);
                gemState_.iterLast = gemState_.iterGlobal;
                return;
            }
        }
    }

    for (int i = 0; i < state_.nConPhases; ++i) {
        if (state_.dMolesPhase(i) < 0 ||
            state_.dMolesPhase(i) < state_.tolerances[kTolPhaseMoles]) {
            int speciesIdx = state_.iAssemblage(i);
            if (speciesIdx > 0) {
                removePureConPhase(speciesIdx - 1);
                gemState_.iterLast = gemState_.iterGlobal;
                return;
            }
        }
    }

    // Check for phases that should be added (positive driving force)
    // Set default step size if not set
    if (gemState_.iterStep == 0) {
        gemState_.iterStep = 10;
    }

    // Check phase stability - allow Newton solver to converge before adding phases
    int checkInterval = (gemState_.iterGlobal < 100) ? 50 : gemState_.iterStep * 5;
    if (gemState_.iterGlobal - gemState_.iterLast < checkInterval) {
        return;
    }

    // Compute driving forces for all unstable solution phases
    computeSolutionPhaseDrivingForces();

    // Find phase with highest driving force (most favorable to add)
    int bestPhase = -1;
    double bestDF = state_.tolerances[kTolDrivingForce];
    for (int iPhase = 0; iPhase < state_.nSolnPhasesSys; ++iPhase) {
        if (gemState_.lSolnPhases[iPhase]) {
            continue;  // Already stable
        }

        double drivingForce = gemState_.dDrivingForceSoln(iPhase);
        if (drivingForce > bestDF) {
            bestDF = drivingForce;
            bestPhase = iPhase;
        }
    }

    // Try to add the best phase
    if (bestPhase >= 0) {
        if (addSolnPhase(bestPhase)) {
            gemState_.iterLast = gemState_.iterGlobal;
            gemState_.iterLastSoln = gemState_.iterGlobal;
            return;
        }

        // If couldn't add due to phase rule, try swapping with worst pure phase
        if (state_.nConPhases > 0) {
            int worstPhaseIdx = -1;
            double worstDF = -1e30;

            for (int i = 0; i < state_.nConPhases; ++i) {
                int speciesIdx = state_.iAssemblage(i) - 1;
                if (speciesIdx < 0 || speciesIdx >= state_.nSpecies) continue;

                double muCalc = 0.0;
                for (int j = 0; j < state_.nElements; ++j) {
                    muCalc += state_.dElementPotential(j) * state_.dStoichSpecies(speciesIdx, j);
                }
                muCalc /= state_.iParticlesPerMole(speciesIdx);

                double phaseDF = state_.dStdGibbsEnergy(speciesIdx) - muCalc;

                if (phaseDF > worstDF) {
                    worstDF = phaseDF;
                    worstPhaseIdx = i;
                }
            }

            if (worstPhaseIdx >= 0 && bestDF > -worstDF) {
                int removedSpecies = state_.iAssemblage(worstPhaseIdx) - 1;
                removePureConPhase(removedSpecies);

                if (addSolnPhase(bestPhase)) {
                    gemState_.iterLast = gemState_.iterGlobal;
                    gemState_.iterLastSoln = gemState_.iterGlobal;
                    return;
                }
            }
        }
    }

    // Check pure condensed phases
    for (int i = state_.nSpeciesPhase(state_.nSolnPhasesSys); i < state_.nSpecies; ++i) {
        double muCalc = 0.0;
        for (int j = 0; j < state_.nElements; ++j) {
            muCalc += state_.dElementPotential(j) * state_.dStoichSpecies(i, j);
        }
        muCalc /= state_.iParticlesPerMole(i);

        double drivingForce = state_.dStdGibbsEnergy(i) - muCalc;

        if (drivingForce < -state_.tolerances[kTolDrivingForce]) {
            if (addPureConPhase(i)) {
                gemState_.iterLast = gemState_.iterGlobal;
                gemState_.iterLastCon = gemState_.iterGlobal;
                return;
            }
        }
    }
}

bool PhaseAssemblageManager::addSolnPhase(int phaseIndex) {
    // Count active elements
    int nActiveElements = 0;
    for (int j = 0; j < state_.nElements - state_.nChargedConstraints; ++j) {
        if (state_.dMolesElement(j) > 0.0) {
            ++nActiveElements;
        }
    }

    // Check Gibbs phase rule: P <= C + 2
    // For isothermal-isobaric: P <= C
    int maxPhases = nActiveElements;
    int currentPhases = state_.nConPhases + state_.nSolnPhases;

    if (currentPhases >= maxPhases) {
        return false;  // Would violate phase rule
    }

    // Check if phase is feasible (all species contain only active elements)
    int iFirst = (phaseIndex > 0) ? state_.nSpeciesPhase(phaseIndex) : 0;
    int iLast = state_.nSpeciesPhase(phaseIndex + 1);

    bool hasFeasibleSpecies = false;
    for (int i = iFirst; i < iLast; ++i) {
        if (speciesIsFeasible(i)) {
            hasFeasibleSpecies = true;
            break;
        }
    }

    if (!hasFeasibleSpecies) {
        return false;
    }

    // Save current assemblage
    saveCurrentAssemblage();

    // Add the phase
    ++state_.nSolnPhases;
    int newPhaseIdx = state_.nElements - state_.nSolnPhases;
    state_.iAssemblage(newPhaseIdx) = -(phaseIndex + 1);
    gemState_.lSolnPhases[phaseIndex] = true;

    // Initialize phase moles
    double initialMoles = 1e-6;
    state_.dMolesPhase(newPhaseIdx) = initialMoles;

    // Count feasible species in this phase
    int nFeasibleSpecies = 0;
    for (int i = iFirst; i < iLast; ++i) {
        if (speciesIsFeasible(i)) {
            ++nFeasibleSpecies;
        }
    }

    // Initialize species mole fractions for this phase (normalize by feasible species only)
    if (nFeasibleSpecies > 0) {
        double fracPerSpecies = 1.0 / nFeasibleSpecies;
        for (int i = iFirst; i < iLast; ++i) {
            if (speciesIsFeasible(i)) {
                state_.dMolFraction(i) = fracPerSpecies;
                state_.dMolesSpecies(i) = fracPerSpecies * initialMoles;
            } else {
                state_.dMolFraction(i) = 0.0;
                state_.dMolesSpecies(i) = 0.0;
            }
        }
    } else {
        // No feasible species - phase should not have been added
        for (int i = iFirst; i < iLast; ++i) {
            state_.dMolFraction(i) = 0.0;
            state_.dMolesSpecies(i) = 0.0;
        }
    }

    return true;
}

bool PhaseAssemblageManager::removeSolnPhase(int phaseIndex) {
    if (phaseIndex < 0 || phaseIndex >= state_.nSolnPhases) {
        return false;
    }

    // Save current assemblage
    saveCurrentAssemblage();

    // Get system phase index
    int systemPhaseIdx = -state_.iAssemblage(state_.nElements - state_.nSolnPhases + phaseIndex) - 1;
    if (systemPhaseIdx < 0) {
        return false;
    }

    // Mark phase as unstable
    gemState_.lSolnPhases[systemPhaseIdx] = false;

    // Zero out species in this phase
    int iFirst = (systemPhaseIdx > 0) ? state_.nSpeciesPhase(systemPhaseIdx) : 0;
    int iLast = state_.nSpeciesPhase(systemPhaseIdx + 1);
    for (int i = iFirst; i < iLast; ++i) {
        state_.dMolesSpecies(i) = 0.0;
        state_.dMolFraction(i) = 0.0;
    }

    // Shift assemblage entries
    for (int i = phaseIndex; i < state_.nSolnPhases - 1; ++i) {
        int idx = state_.nElements - state_.nSolnPhases + i;
        state_.iAssemblage(idx) = state_.iAssemblage(idx + 1);
        state_.dMolesPhase(idx) = state_.dMolesPhase(idx + 1);
    }

    // Clear the last solution phase slot to prevent stale data
    int lastSlot = state_.nElements - 1;
    state_.iAssemblage(lastSlot) = 0;
    state_.dMolesPhase(lastSlot) = 0.0;

    --state_.nSolnPhases;
    return true;
}

bool PhaseAssemblageManager::addPureConPhase(int speciesIndex) {
    // Check if species is feasible
    if (!speciesIsFeasible(speciesIndex)) {
        return false;
    }

    // Prevent duplicate pure condensed phases
    for (int i = 0; i < state_.nConPhases; ++i) {
        if (state_.iAssemblage(i) == speciesIndex + 1) {
            return false;
        }
    }

    // Count active elements
    int nActiveElements = 0;
    for (int j = 0; j < state_.nElements - state_.nChargedConstraints; ++j) {
        if (state_.dMolesElement(j) > 0.0) {
            ++nActiveElements;
        }
    }

    // Check phase rule
    int maxPhases = nActiveElements;
    int currentPhases = state_.nConPhases + state_.nSolnPhases;

    if (currentPhases >= maxPhases) {
        return false;
    }

    // Save current assemblage
    saveCurrentAssemblage();

    // Add the phase
    ++state_.nConPhases;
    double initialMoles = 1e-6;
    state_.iAssemblage(state_.nConPhases - 1) = speciesIndex + 1;
    state_.dMolesPhase(state_.nConPhases - 1) = initialMoles;

    // Set species moles and mole fraction consistently for pure condensed phase
    state_.dMolesSpecies(speciesIndex) = initialMoles;
    state_.dMolFraction(speciesIndex) = 1.0;

    return true;
}

bool PhaseAssemblageManager::removePureConPhase(int speciesIndex) {
    // Find the phase in the assemblage
    int phaseIdx = -1;
    for (int i = 0; i < state_.nConPhases; ++i) {
        if (state_.iAssemblage(i) == speciesIndex + 1) {
            phaseIdx = i;
            break;
        }
    }

    if (phaseIdx < 0) {
        return false;
    }

    // Save current assemblage
    saveCurrentAssemblage();

    // Clear species moles for the removed pure phase
    state_.dMolesSpecies(speciesIndex) = 0.0;
    state_.dMolFraction(speciesIndex) = 0.0;

    // Shift assemblage entries
    for (int i = phaseIdx; i < state_.nConPhases - 1; ++i) {
        state_.iAssemblage(i) = state_.iAssemblage(i + 1);
        state_.dMolesPhase(i) = state_.dMolesPhase(i + 1);
    }

    // Clear the last condensed phase slot to prevent stale data
    state_.iAssemblage(state_.nConPhases - 1) = 0;
    state_.dMolesPhase(state_.nConPhases - 1) = 0.0;

    --state_.nConPhases;
    return true;
}

void PhaseAssemblageManager::saveCurrentAssemblage() {
    // Save phase moles for potential revert
    gemState_.dMolesPhaseLast = state_.dMolesPhase;

    // Update iteration history with current assemblage
    // This is critical for revert() to restore the correct assemblage
    for (int i = 0; i < state_.nElements; ++i) {
        gemState_.iterHistory(i, gemState_.iterGlobal) = state_.iAssemblage(i);
    }
}

void PhaseAssemblageManager::revert() {
    // Revert to a previous successful assemblage using iteration history
    int revertIter = std::max(1, gemState_.iterGlobal - 200);

    if (gemState_.iterGlobal > 1500) {
        revertIter = std::max(1, gemState_.iterGlobal - 1000);
    }

    // Restore assemblage from history
    for (int i = 0; i < state_.nElements; ++i) {
        state_.iAssemblage(i) = gemState_.iterHistory(i, revertIter);
    }

    // Recalculate phase counts
    state_.nConPhases = 0;
    state_.nSolnPhases = 0;

    for (int i = 0; i < state_.nElements; ++i) {
        int idx = state_.iAssemblage(i);
        if (idx > 0) {
            state_.nConPhases++;
        } else if (idx < 0) {
            state_.nSolnPhases++;
        }
    }

    // Update stable phase flags
    std::fill(gemState_.lSolnPhases.begin(), gemState_.lSolnPhases.end(), false);
    for (int i = 0; i < state_.nElements; ++i) {
        int idx = state_.iAssemblage(i);
        if (idx < 0) {
            int phaseIdx = -idx - 1;
            if (phaseIdx >= 0 && phaseIdx < static_cast<int>(gemState_.lSolnPhases.size())) {
                gemState_.lSolnPhases[phaseIdx] = true;
            }
        }
    }

    gemState_.iterRevert = gemState_.iterGlobal;
    gemState_.dGEMFunctionNorm *= 10.0;  // Relax convergence tolerance
}

} // namespace Thermochimica
