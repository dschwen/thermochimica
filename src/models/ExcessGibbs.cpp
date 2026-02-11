#include "thermochimica/ThermoContext.hpp"
#include "thermochimica/models/ModelFactory.hpp"

namespace Thermochimica {

// Forward declarations for legacy excess Gibbs energy models
void computeExcessGibbsQKTO(ThermoContext& ctx, int phaseIndex);
void computeExcessGibbsRKMP(ThermoContext& ctx, int phaseIndex);
void computeExcessGibbsSUBL(ThermoContext& ctx, int phaseIndex);
void computeExcessGibbsSUBG(ThermoContext& ctx, int phaseIndex);

// Global model factory instance (initialized on first use)
static ModelFactory& getModelFactory() {
    static ModelFactory factory;
    return factory;
}

// NEW: Factory-based excess Gibbs computation using IThermodynamicModel
// This is the preferred method going forward
void computeExcessGibbsWithFactory(ThermoContext& ctx, int phaseIndex) {
    auto& thermo = *ctx.thermo;
    auto& gem = *ctx.gem;

    if (phaseIndex < 0 || phaseIndex >= static_cast<int>(thermo.iSolnPhaseType.size())) {
        return;
    }

    Constants::PhaseType phaseType = thermo.iSolnPhaseType[phaseIndex];

    // Get model from factory
    IThermodynamicModel* model = getModelFactory().getModel(phaseType);

    if (model != nullptr && model->canHandle(thermo, phaseIndex)) {
        // Use factory model
        model->computeExcessGibbs(thermo, gem, phaseIndex);
    } else {
        // Fallback: unimplemented model (set excess Gibbs to zero)
        thermo.dGibbsSolnPhase(phaseIndex) = 0.0;
    }
}

// LEGACY: Dispatcher-based excess Gibbs computation (backward compatibility)
// Routes to appropriate legacy model based on phase type
void computeExcessGibbs(ThermoContext& ctx, int phaseIndex) {
    // Use new factory-based implementation
    computeExcessGibbsWithFactory(ctx, phaseIndex);

    /* LEGACY DISPATCHER (commented out, replaced by factory)
    auto& thermo = *ctx.thermo;

    if (phaseIndex < 0 || phaseIndex >= static_cast<int>(thermo.iSolnPhaseType.size())) {
        return;
    }

    Constants::PhaseType phaseType = thermo.iSolnPhaseType[phaseIndex];

    if (phaseType == Constants::PhaseType::IDMX) {
        // Ideal mixing - no excess Gibbs energy
        thermo.dGibbsSolnPhase(phaseIndex) = 0.0;
    } else if (phaseType == Constants::PhaseType::QKTO) {
        // Kohler-Toop model
        // Note: QKTO model uses 1-based phase indexing internally
        computeExcessGibbsQKTO(ctx, phaseIndex + 1);
    } else if (phaseType == Constants::PhaseType::RKMP || phaseType == Constants::PhaseType::RKMPM) {
        // Redlich-Kister-Muggianu
        // computeExcessGibbsRKMP(ctx, phaseIndex);
    } else if (phaseType == Constants::PhaseType::SUBL || phaseType == Constants::PhaseType::SUBLM) {
        // Compound Energy Formalism
        // computeExcessGibbsSUBL(ctx, phaseIndex);
    } else if (phaseType == Constants::PhaseType::SUBG || phaseType == Constants::PhaseType::SUBQ) {
        // Modified Quasichemical Model
        // computeExcessGibbsSUBG(ctx, phaseIndex);
    } else if (phaseType == Constants::PhaseType::SUBI) {
        // Ionic sublattice
        // computeExcessGibbsSUBI(ctx, phaseIndex);
    }
    */
}

// Wrapper function for compatibility with Subminimization
// This is the function signature expected by the Fortran-style code
void compExcessGibbsEnergy(ThermoContext& ctx, int phaseIndex) {
    // NOTE: Subminimization passes 1-based index, but computeExcessGibbs expects 0-based
    // So we subtract 1 here to convert
    computeExcessGibbs(ctx, phaseIndex - 1);
}

} // namespace Thermochimica
