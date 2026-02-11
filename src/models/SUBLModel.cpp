/// @file SUBLModel.cpp
/// @brief Implementation of SUBL thermodynamic model

#include "thermochimica/models/SUBLModel.hpp"
#include "thermochimica/context/ThermoState.hpp"
#include "thermochimica/context/GEMState.hpp"
#include "thermochimica/ThermoContext.hpp"

namespace Thermochimica {

// Forward declaration of legacy SUBL function
void compExcessGibbsEnergySUBL(ThermoContext& ctx, int iSolnIndex);

void SUBLModel::computeExcessGibbs(ThermoState& state, GEMState& gemState, int phaseIndex) {
    // Create temporary context for bridge to legacy implementation
    ThermoContext ctx;

    // Move state references into context (non-owning)
    ThermoState* statePtr = &state;
    GEMState* gemPtr = &gemState;

    ctx.thermo.reset(statePtr);
    ctx.gem.reset(gemPtr);

    // Call legacy SUBL implementation (uses 1-based phase indexing)
    compExcessGibbsEnergySUBL(ctx, phaseIndex + 1);

    // Release pointers so they won't be deleted when ctx goes out of scope
    ctx.thermo.release();
    ctx.gem.release();
}

bool SUBLModel::canHandle(const ThermoState& state, int phaseIndex) const {
    if (phaseIndex < 0 || phaseIndex >= static_cast<int>(state.iSolnPhaseType.size())) {
        return false;
    }

    Constants::PhaseType phaseType = state.iSolnPhaseType[phaseIndex];
    return (phaseType == Constants::PhaseType::SUBL ||
            phaseType == Constants::PhaseType::SUBLM);
}

} // namespace Thermochimica
