#include "thermochimica/solver/GEMSolver.hpp"
#include "thermochimica/util/Constants.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

namespace Thermochimica {

void GEMLineSearch::search(ThermoContext& ctx) {
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;

    double stepLength = 1.0;

    // Save current state before line search
    Eigen::VectorXd savedElementPotential = thermo.dElementPotential;
    Eigen::VectorXd savedMolesPhase = thermo.dMolesPhase;
    Eigen::VectorXd savedMolesSpecies = thermo.dMolesSpecies;

    // Initialize line search
    init(ctx, stepLength);

    // Save current function norm
    double functionNormLast = gem.dGEMFunctionNorm;

    // Wolfe line search loop (max 5 iterations)
    for (int iter = 0; iter < 5; ++iter) {
        // Restore saved state before each trial
        thermo.dElementPotential = savedElementPotential;
        thermo.dMolesPhase = savedMolesPhase;
        thermo.dMolesSpecies = savedMolesSpecies;

        // Update element potentials
        updateElementPotentials(ctx, stepLength);

        // Update phase moles
        updatePhaseMoles(ctx, stepLength);

        // Update species moles and fractions
        updateSpecies(ctx);

        // Compute new function norm (mass balance residual)
        double newNorm = 0.0;
        for (int j = 0; j < thermo.nElements; ++j) {
            double residual = thermo.dMolesElement(j);

            // Subtract species contributions
            for (int i = 0; i < thermo.nSpecies; ++i) {
                if (thermo.dMolesSpecies(i) > 0) {
                    residual -= thermo.dMolesSpecies(i) * thermo.dStoichSpecies(i, j);
                }
            }

            newNorm += residual * residual;
        }
        gem.dGEMFunctionNorm = std::sqrt(newNorm);

        // Check Wolfe conditions
        if (checkWolfeConditions(ctx, functionNormLast)) {
            break;
        }

        // Dampen step length
        stepLength *= 0.5;
    }

    // Update last values
    gem.dGEMFunctionNormLast = gem.dGEMFunctionNorm;
    gem.dMolesPhaseLast = thermo.dMolesPhase;
}

void GEMLineSearch::init(ThermoContext& ctx, double& stepLength) {
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;

    stepLength = 1.0;

    // Constrain step length based on iteration count
    double maxGamma = 1.0;
    if (gem.iterGlobal < 1000) {
        maxGamma = 1.0;
    } else if (gem.iterGlobal < 2000) {
        maxGamma = 0.5;
    } else {
        maxGamma = 0.25;
    }

    // dUpdateVar is the Newton step (delta), not target values
    // Limit step based on element potential changes
    for (int i = 0; i < thermo.nElements; ++i) {
        double change = std::abs(gem.dUpdateVar(i));  // This is the delta
        if (change > 1e-10) {
            double limit = maxGamma / change;
            stepLength = std::min(stepLength, limit);
        }
    }

    // Limit step based on phase mole changes
    int nElements = thermo.nElements;
    int nSolnPhases = thermo.nSolnPhases;

    for (int i = 0; i < nSolnPhases; ++i) {
        int idx = nElements + i;
        if (idx < gem.dUpdateVar.size()) {
            double molesCurrent = thermo.dMolesPhase(nElements - nSolnPhases + i);
            double molesDelta = gem.dUpdateVar(idx);  // This is the delta

            // Prevent phase moles from going negative
            // new = current + step * delta
            // want: current + step * delta >= 0
            // if delta < 0: step <= current / |delta|
            if (molesDelta < 0 && molesCurrent > 0) {
                double limit = 0.9 * molesCurrent / std::abs(molesDelta);
                stepLength = std::min(stepLength, limit);
            }

            // Limit large changes (both increase and decrease)
            if (molesCurrent > 0) {
                double maxChange = 0.3 * molesCurrent;  // Allow at most 30% change
                if (std::abs(molesDelta) > maxChange) {
                    double limit = maxChange / std::abs(molesDelta);
                    stepLength = std::min(stepLength, limit);
                }
            }

            // Also limit absolute growth to prevent runaway
            // Don't allow phase moles to more than double in one iteration
            if (molesDelta > 0) {
                double maxGrowth = std::max(molesCurrent, 0.1);  // At most double, min 0.1
                double limit = maxGrowth / molesDelta;
                stepLength = std::min(stepLength, limit);
            }
        }
    }

    // Absolute maximum step length constraint
    stepLength = std::min(stepLength, 0.1);

    stepLength = std::max(stepLength, 1e-10);
}

bool GEMLineSearch::checkWolfeConditions(ThermoContext& ctx, double functionNormLast) {
    auto& gem = *ctx.gem;

    // Sufficient decrease condition
    if (gem.dGEMFunctionNorm < 1e-6) {
        return true;  // Excellent progress
    }

    // Check relative change
    double relChange = gem.dGEMFunctionNorm / std::max(functionNormLast, 1e-20);

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

void GEMLineSearch::updateElementPotentials(ThermoContext& ctx, double stepLength) {
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;

    // dUpdateVar contains the Newton step (delta), not the target value
    // Apply: new = old + stepLength * delta
    for (int i = 0; i < thermo.nElements; ++i) {
        double oldVal = thermo.dElementPotential(i);
        double newVal = oldVal + stepLength * gem.dUpdateVar(i);
        thermo.dElementPotential(i) = newVal;
    }
}

void GEMLineSearch::updatePhaseMoles(ThermoContext& ctx, double stepLength) {
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;

    int nElements = thermo.nElements;
    int nSolnPhases = thermo.nSolnPhases;
    int nConPhases = thermo.nConPhases;

    // dUpdateVar contains the Newton step (delta), not the target value
    // Apply: new = old + stepLength * delta

    // Update solution phase moles
    for (int i = 0; i < nSolnPhases; ++i) {
        int updateIdx = nElements + i;
        int phaseIdx = nElements - nSolnPhases + i;

        if (updateIdx < gem.dUpdateVar.size() && phaseIdx < thermo.dMolesPhase.size()) {
            double oldMoles = thermo.dMolesPhase(phaseIdx);
            double newMoles = oldMoles + stepLength * gem.dUpdateVar(updateIdx);
            thermo.dMolesPhase(phaseIdx) = std::max(newMoles, 0.0);
        }
    }

    // Update pure condensed phase moles
    for (int i = 0; i < nConPhases; ++i) {
        int updateIdx = nElements + nSolnPhases + i;

        if (updateIdx < gem.dUpdateVar.size() && i < thermo.dMolesPhase.size()) {
            double oldMoles = thermo.dMolesPhase(i);
            double newMoles = oldMoles + stepLength * gem.dUpdateVar(updateIdx);
            thermo.dMolesPhase(i) = std::max(newMoles, 0.0);
        }
    }
}

void GEMLineSearch::updateSpecies(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;

    // Keep mole fractions FIXED for now.
    // The mole fractions were initialized to match the overall system composition.
    // Only update species moles = phase moles * mole fraction.
    //
    // This is a simplification - a proper GEM solver would update mole fractions
    // based on element potentials, but that requires careful handling to maintain
    // mass balance.

    for (int iPhase = 0; iPhase < thermo.nSolnPhases; ++iPhase) {
        int assembIdx = thermo.nElements - thermo.nSolnPhases + iPhase;
        int phaseIdx = -thermo.iAssemblage(assembIdx) - 1;
        if (phaseIdx < 0 || phaseIdx >= thermo.nSolnPhasesSys) continue;

        int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
        int iLast = thermo.nSpeciesPhase(phaseIdx + 1);
        double phaseMoles = thermo.dMolesPhase(assembIdx);

        // Update species moles from phase moles and FIXED mole fractions
        for (int i = iFirst; i < iLast; ++i) {
            if (thermo.dMolFraction(i) > 0) {
                thermo.dMolesSpecies(i) = phaseMoles * thermo.dMolFraction(i);
            }
        }
    }
}

} // namespace Thermochimica
