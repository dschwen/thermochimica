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

    // Temporarily move state references into context
    ThermoState* statePtr = &state;
    ThermoIO* ioPtr = &io;
    GEMState* gemPtr = &gemState;
    PhaseConstraints* constraintsPtr = &constraints_;

    ctx.thermo.reset(statePtr);
    ctx.io.reset(ioPtr);
    ctx.gem.reset(gemPtr);
    ctx.phaseConstraints.reset(constraintsPtr);

    // Call legacy solver (will detect active constraints and run constrained path)
    int result = GEMSolver::solve(ctx);

    // Release pointers so they won't be deleted when ctx goes out of scope
    ctx.thermo.release();
    ctx.io.release();
    ctx.gem.release();
    ctx.phaseConstraints.release();

    return result;
}

void ConstrainedGEMSolver::initialize(ThermoState& state, GEMState& gemState) {
    // Create temporary context for bridge to legacy init
    ThermoContext ctx;

    ThermoState* statePtr = &state;
    GEMState* gemPtr = &gemState;

    ctx.thermo.reset(statePtr);
    ctx.gem.reset(gemPtr);

    // Call legacy init
    GEMSolver::init(ctx);

    // Release pointers
    ctx.thermo.release();
    ctx.gem.release();
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
