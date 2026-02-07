/// @file WolfeLineSearch.cpp
/// @brief Implementation of Wolfe conditions line search

#include "thermochimica/solver/WolfeLineSearch.hpp"
#include "thermochimica/context/ThermoState.hpp"
#include "thermochimica/context/GEMState.hpp"
#include <cmath>
#include <algorithm>

namespace Thermochimica {

void WolfeLineSearch::search(ThermoState& state,
                             GEMState& gemState,
                             const Eigen::VectorXd& direction,
                             double& stepLength) {
    // Save current state before line search
    Eigen::VectorXd savedElementPotential = state.dElementPotential;
    Eigen::VectorXd savedMolesPhase = state.dMolesPhase;
    Eigen::VectorXd savedMolesSpecies = state.dMolesSpecies;

    // Store direction in gemState for helper methods (temporary bridge)
    gemState.dUpdateVar = direction;

    // Initialize step length with constraints
    initStepLength(state, gemState, stepLength);

    // Save current function norm
    double functionNormLast = gemState.dGEMFunctionNorm;

    // Wolfe line search loop (max 5 iterations)
    for (int iter = 0; iter < 5; ++iter) {
        // Restore saved state before each trial
        state.dElementPotential = savedElementPotential;
        state.dMolesPhase = savedMolesPhase;
        state.dMolesSpecies = savedMolesSpecies;

        // Update state with current step length
        updateState(state, gemState, direction, stepLength);

        // Compute new function norm (mass balance residual)
        double newNorm = 0.0;
        for (int j = 0; j < state.nElements; ++j) {
            double residual = state.dMolesElement(j);

            // Subtract species contributions
            for (int i = 0; i < state.nSpecies; ++i) {
                if (state.dMolesSpecies(i) > 0) {
                    residual -= state.dMolesSpecies(i) * state.dStoichSpecies(i, j);
                }
            }

            newNorm += residual * residual;
        }
        gemState.dGEMFunctionNorm = std::sqrt(newNorm);

        // Check Wolfe conditions
        if (checkWolfeConditions(gemState, functionNormLast)) {
            break;
        }

        // Dampen step length
        stepLength *= 0.5;
    }

    // Update last values
    gemState.dGEMFunctionNormLast = gemState.dGEMFunctionNorm;
    gemState.dMolesPhaseLast = state.dMolesPhase;
}

void WolfeLineSearch::updateState(ThermoState& state,
                                  GEMState& gemState,
                                  const Eigen::VectorXd& direction,
                                  double stepLength) {
    // Update element potentials
    updateElementPotentials(state, gemState, direction, stepLength);

    // Update phase moles
    updatePhaseMoles(state, gemState, direction, stepLength);

    // Update species moles and fractions
    updateSpecies(state);
}

void WolfeLineSearch::initStepLength(ThermoState& state,
                                     GEMState& gemState,
                                     double& stepLength) {
    stepLength = 1.0;

    // Constrain step length based on iteration count
    double maxGamma = 1.0;
    if (gemState.iterGlobal < 1000) {
        maxGamma = 1.0;
    } else if (gemState.iterGlobal < 2000) {
        maxGamma = 0.5;
    } else {
        maxGamma = 0.25;
    }

    const Eigen::VectorXd& direction = gemState.dUpdateVar;

    // Limit step based on element potential changes
    for (int i = 0; i < state.nElements; ++i) {
        double change = std::abs(direction(i));
        if (change > 1e-10) {
            double limit = maxGamma / change;
            stepLength = std::min(stepLength, limit);
        }
    }

    // Limit step based on phase mole changes
    int nElements = state.nElements;
    int nSolnPhases = state.nSolnPhases;

    for (int i = 0; i < nSolnPhases; ++i) {
        int idx = nElements + i;
        if (idx < direction.size()) {
            double molesCurrent = state.dMolesPhase(nElements - nSolnPhases + i);
            double molesDelta = direction(idx);

            // Prevent phase moles from going negative
            if (molesDelta < 0 && molesCurrent > 0) {
                double limit = 0.9 * molesCurrent / std::abs(molesDelta);
                stepLength = std::min(stepLength, limit);
            }

            // Limit large changes (at most 30% change)
            if (molesCurrent > 0) {
                double maxChange = 0.3 * molesCurrent;
                if (std::abs(molesDelta) > maxChange) {
                    double limit = maxChange / std::abs(molesDelta);
                    stepLength = std::min(stepLength, limit);
                }
            }

            // Limit absolute growth (don't allow phase to more than double)
            if (molesDelta > 0) {
                double maxGrowth = std::max(molesCurrent, 0.1);
                double limit = maxGrowth / molesDelta;
                stepLength = std::min(stepLength, limit);
            }
        }
    }

    // Absolute maximum step length constraint
    stepLength = std::min(stepLength, 0.1);
    stepLength = std::max(stepLength, 1e-10);
}

bool WolfeLineSearch::checkWolfeConditions(const GEMState& gemState,
                                           double functionNormLast) const {
    // Sufficient decrease condition
    if (gemState.dGEMFunctionNorm < 1e-6) {
        return true;  // Excellent progress
    }

    // Check relative change
    double relChange = gemState.dGEMFunctionNorm / std::max(functionNormLast, 1e-20);

    if (relChange < 0.999) {
        return true;  // Sufficient decrease
    }

    if (relChange > 0.97 && relChange < 1.0) {
        return true;  // Plateau detected
    }

    if (relChange >= 1.0) {
        return true;  // Diverging, exit and let phase assemblage adjust
    }

    return false;
}

void WolfeLineSearch::updateElementPotentials(ThermoState& state,
                                              GEMState& gemState,
                                              const Eigen::VectorXd& direction,
                                              double stepLength) {
    // Apply: new = old + stepLength * delta
    for (int i = 0; i < state.nElements; ++i) {
        double oldVal = state.dElementPotential(i);
        double newVal = oldVal + stepLength * direction(i);
        state.dElementPotential(i) = newVal;
    }
}

void WolfeLineSearch::updatePhaseMoles(ThermoState& state,
                                       GEMState& gemState,
                                       const Eigen::VectorXd& direction,
                                       double stepLength) {
    int nElements = state.nElements;
    int nSolnPhases = state.nSolnPhases;
    int nConPhases = state.nConPhases;

    // Update solution phase moles
    for (int i = 0; i < nSolnPhases; ++i) {
        int updateIdx = nElements + i;
        int phaseIdx = nElements - nSolnPhases + i;

        if (updateIdx < direction.size() && phaseIdx < state.dMolesPhase.size()) {
            double oldMoles = state.dMolesPhase(phaseIdx);
            double newMoles = oldMoles + stepLength * direction(updateIdx);
            state.dMolesPhase(phaseIdx) = std::max(newMoles, 0.0);
        }
    }

    // Update pure condensed phase moles
    for (int i = 0; i < nConPhases; ++i) {
        int updateIdx = nElements + nSolnPhases + i;

        if (updateIdx < direction.size() && i < state.dMolesPhase.size()) {
            double oldMoles = state.dMolesPhase(i);
            double newMoles = oldMoles + stepLength * direction(updateIdx);
            state.dMolesPhase(i) = std::max(newMoles, 0.0);
        }
    }
}

void WolfeLineSearch::updateSpecies(ThermoState& state) {
    // Update species moles from phase moles and mole fractions
    // Note: This keeps mole fractions FIXED during line search

    for (int iPhase = 0; iPhase < state.nSolnPhases; ++iPhase) {
        int assembIdx = state.nElements - state.nSolnPhases + iPhase;
        int phaseIdx = -state.iAssemblage(assembIdx) - 1;
        if (phaseIdx < 0 || phaseIdx >= state.nSolnPhasesSys) continue;

        int iFirst = (phaseIdx > 0) ? state.nSpeciesPhase(phaseIdx) : 0;
        int iLast = state.nSpeciesPhase(phaseIdx + 1);
        double phaseMoles = state.dMolesPhase(assembIdx);

        // Update species moles from phase moles and FIXED mole fractions
        for (int i = iFirst; i < iLast; ++i) {
            if (state.dMolFraction(i) > 0) {
                state.dMolesSpecies(i) = phaseMoles * state.dMolFraction(i);
            }
        }
    }
}

} // namespace Thermochimica
