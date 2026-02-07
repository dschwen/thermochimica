/// @file StandardGEMSolver.cpp
/// @brief Implementation of standard GEM solver

#include "thermochimica/solver/StandardGEMSolver.hpp"
#include "thermochimica/solver/GEMSolver.hpp"
#include "thermochimica/context/ThermoState.hpp"
#include "thermochimica/context/ThermoIO.hpp"
#include "thermochimica/context/GEMState.hpp"
#include "thermochimica/context/PhaseConstraints.hpp"
#include "thermochimica/ThermoContext.hpp"

namespace Thermochimica {

StandardGEMSolver::StandardGEMSolver()
    : phaseConstraints_(std::make_unique<PhaseConstraints>()) {
}

int StandardGEMSolver::solve(ThermoState& state,
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

    ctx.thermo.reset(statePtr);
    ctx.io.reset(ioPtr);
    ctx.gem.reset(gemPtr);

    // Transfer PhaseConstraints ownership to context for unconstrained solve
    // (will be transferred back after solve completes)
    ctx.phaseConstraints.reset(phaseConstraints_.release());

    // Call legacy solver
    int result = GEMSolver::solve(ctx);

    // Release pointers so they won't be deleted when ctx goes out of scope
    ctx.thermo.release();
    ctx.io.release();
    ctx.gem.release();

    // Transfer PhaseConstraints ownership back to StandardGEMSolver
    phaseConstraints_.reset(ctx.phaseConstraints.release());

    return result;
}

void StandardGEMSolver::initialize(ThermoState& state, GEMState& gemState) {
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

bool StandardGEMSolver::isConverged(const ThermoState& state,
                                   const GEMState& gemState) const {
    return gemState.lConverged;
}

} // namespace Thermochimica
