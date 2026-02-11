/// @file ConstrainedGEMSolver.cpp
/// @brief Implementation of constrained GEM solver

#include "thermochimica/solver/ConstrainedGEMSolver.hpp"
#include "thermochimica/solver/GEMSolver.hpp"
#include "thermochimica/context/ThermoState.hpp"
#include "thermochimica/context/ThermoIO.hpp"
#include "thermochimica/context/GEMState.hpp"
#include "thermochimica/context/PhaseConstraints.hpp"
#include "thermochimica/ThermoContext.hpp"

namespace Thermochimica {

ConstrainedGEMSolver::ConstrainedGEMSolver(PhaseConstraints& constraints)
    : constraints_(constraints) {
}

int ConstrainedGEMSolver::solve(ThermoState& state,
                                ThermoIO& io,
                                GEMState& gemState,
                                PhaseAssemblageManager& phaseManager,
                                INewtonSolver& newton,
                                ILineSearch& lineSearch,
                                const std::vector<IThermodynamicModel*>& models) {
    // Create temporary context for bridge to legacy solver
    ThermoContext ctx;

    // Temporarily borrow state references
    ctx.thermo.reset(&state);
    ctx.io.reset(&io);
    ctx.gem.reset(&gemState);
    ctx.phaseConstraints.reset(&constraints_);

    // RAII guard ensures release() is called even if exception is thrown
    struct Guard {
        ThermoContext& ctx;
        ~Guard() {
            ctx.thermo.release();
            ctx.io.release();
            ctx.gem.release();
            ctx.phaseConstraints.release();
        }
    } guard{ctx};

    // Call legacy solver (exception-safe with guard)
    int result = GEMSolver::solve(ctx);

    return result;
}

void ConstrainedGEMSolver::initialize(ThermoState& state, GEMState& gemState) {
    // Create temporary context for bridge to legacy init
    ThermoContext ctx;

    // Temporarily borrow state references
    ctx.thermo.reset(&state);
    ctx.gem.reset(&gemState);

    // RAII guard ensures release() is called even if exception is thrown
    struct Guard {
        ThermoContext& ctx;
        ~Guard() {
            ctx.thermo.release();
            ctx.gem.release();
        }
    } guard{ctx};

    // Call legacy init (exception-safe with guard)
    GEMSolver::init(ctx);
}

bool ConstrainedGEMSolver::isConverged(const ThermoState& state,
                                       const GEMState& gemState) const {
    // For constrained solver, convergence requires both:
    // 1. Inner GEM loop convergence (gemState.lConverged)
    // 2. Phase constraints satisfied
    // The legacy solver sets gemState.lConverged = true only when both are satisfied
    return gemState.lConverged;
}

} // namespace Thermochimica
