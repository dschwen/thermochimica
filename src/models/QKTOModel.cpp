/// @file QKTOModel.cpp
/// @brief Implementation of QKTO thermodynamic model

#include "thermochimica/models/QKTOModel.hpp"
#include "thermochimica/context/ThermoState.hpp"
#include "thermochimica/context/GEMState.hpp"
#include "thermochimica/ThermoContext.hpp"

namespace Thermochimica {

// Forward declaration of existing implementation
void compExcessGibbsEnergyQKTO(ThermoContext& ctx, int phaseIndex);

void QKTOModel::computeExcessGibbs(ThermoState& state,
                                    GEMState& gemState,
                                    int phaseIndex) {
    // Temporary bridge to legacy implementation
    // TODO: Refactor compExcessGibbsEnergyQKTO to work directly with state/gemState

    // Create a temporary context with borrowed pointers
    // We'll manually manage the lifetime to avoid ownership issues
    ThermoContext ctx;

    // Temporarily move the unique_ptrs to avoid deletion
    ThermoState* statePtr = &state;
    GEMState* gemPtr = &gemState;

    ctx.thermo.reset(statePtr);
    ctx.gem.reset(gemPtr);

    // Note: QKTO expects 1-based phase indexing
    compExcessGibbsEnergyQKTO(ctx, phaseIndex + 1);

    // Release the pointers so they won't be deleted when ctx goes out of scope
    ctx.thermo.release();
    ctx.gem.release();
}

bool QKTOModel::canHandle(const ThermoState& thermo, int phaseIndex) const {
    if (phaseIndex < 0 || phaseIndex >= static_cast<int>(thermo.iSolnPhaseType.size())) {
        return false;
    }
    return thermo.iSolnPhaseType[phaseIndex] == Constants::PhaseType::QKTO;
}

} // namespace Thermochimica
