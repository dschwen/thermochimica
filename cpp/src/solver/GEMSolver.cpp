#include "thermochimica/solver/GEMSolver.hpp"
#include "thermochimica/util/ErrorCodes.hpp"
#include "thermochimica/util/Constants.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace Thermochimica {

int GEMSolver::solve(ThermoContext& ctx) {
    auto& io = *ctx.io;
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;

    // Initialize the GEM solver
    init(ctx);

    // Check if system with only pure condensed phases is already converged
    if (thermo.nSolnPhases == 0 && io.INFOThermo == 0 && !io.lReinitLoaded) {
        checkSysOnlyPureConPhases(ctx);
    }

    // Main iteration loop
    for (gem.iterGlobal = 1; gem.iterGlobal <= Constants::kIterGlobalMax; ++gem.iterGlobal) {
        // For systems with solution phases or where mass balance needs adjustment,
        // compute Newton direction and perform line search
        if (thermo.nSolnPhases > 0 || gem.dGEMFunctionNorm > thermo.tolerances[kTolFunctionNorm]) {
            GEMNewton::compute(ctx);
            GEMLineSearch::search(ctx);
        }

        // Check phase assemblage
        PhaseAssemblage::check(ctx);

        // Check convergence
        gem.lConverged = ConvergenceChecker::check(ctx);

        // Exit if error or converged
        if (io.INFOThermo != 0 || gem.lConverged) {
            break;
        }
    }

    // Report error if not converged
    if (!gem.lConverged && io.INFOThermo == 0) {
        io.INFOThermo = ErrorCode::kGEMSolverDidNotConverge;
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
}

} // namespace Thermochimica
