/// @file SUBGModel.cpp
/// @brief Implementation of SUBG thermodynamic model

#include "thermochimica/models/SUBGModel.hpp"
#include "thermochimica/context/ThermoState.hpp"
#include "thermochimica/context/GEMState.hpp"
#include "thermochimica/ThermoContext.hpp"

namespace Thermochimica {

// Forward declaration of legacy SUBG function
void compExcessGibbsEnergySUBG(ThermoContext& ctx, int iSolnIndex);

void SUBGModel::computeExcessGibbs(ThermoState& state, GEMState& gemState, int phaseIndex) {
    // Create temporary context for bridge to legacy implementation
    ThermoContext ctx;

    // Move state references into context (non-owning)
    ThermoState* statePtr = &state;
    GEMState* gemPtr = &gemState;

    ctx.thermo.reset(statePtr);
    ctx.gem.reset(gemPtr);

    // Call legacy SUBG implementation (uses 1-based phase indexing)
    compExcessGibbsEnergySUBG(ctx, phaseIndex + 1);

    // Release pointers so they won't be deleted when ctx goes out of scope
    ctx.thermo.release();
    ctx.gem.release();
}

bool SUBGModel::canHandle(const ThermoState& state, int phaseIndex) const {
    if (phaseIndex < 0 || phaseIndex >= static_cast<int>(state.iSolnPhaseType.size())) {
        return false;
    }

    Constants::PhaseType phaseType = state.iSolnPhaseType[phaseIndex];
    return (phaseType == Constants::PhaseType::SUBG);
}

} // namespace Thermochimica
