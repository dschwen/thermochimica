#include "thermochimica/ThermoContext.hpp"

namespace Thermochimica {

// Forward declarations for excess Gibbs energy models
void computeExcessGibbsQKTO(ThermoContext& ctx, int phaseIndex);
void computeExcessGibbsRKMP(ThermoContext& ctx, int phaseIndex);
void computeExcessGibbsSUBL(ThermoContext& ctx, int phaseIndex);
void computeExcessGibbsSUBG(ThermoContext& ctx, int phaseIndex);

// Excess Gibbs energy computation dispatcher
// Routes to appropriate model based on phase type

void computeExcessGibbs(ThermoContext& ctx, int phaseIndex) {
    auto& thermo = *ctx.thermo;

    if (phaseIndex < 0 || phaseIndex >= static_cast<int>(thermo.cSolnPhaseType.size())) {
        return;
    }

    const std::string& phaseType = thermo.cSolnPhaseType[phaseIndex];

    if (phaseType == "IDMX") {
        // Ideal mixing - no excess Gibbs energy
        thermo.dGibbsSolnPhase(phaseIndex) = 0.0;
    } else if (phaseType == "QKTO") {
        // Kohler-Toop model
        // Note: QKTO model uses 1-based phase indexing internally
        computeExcessGibbsQKTO(ctx, phaseIndex + 1);
    } else if (phaseType == "RKMP" || phaseType == "RKMPM") {
        // Redlich-Kister-Muggianu
        // computeExcessGibbsRKMP(ctx, phaseIndex);
    } else if (phaseType == "SUBL" || phaseType == "SUBLM") {
        // Compound Energy Formalism
        // computeExcessGibbsSUBL(ctx, phaseIndex);
    } else if (phaseType == "SUBG" || phaseType == "SUBQ") {
        // Modified Quasichemical Model
        // computeExcessGibbsSUBG(ctx, phaseIndex);
    } else if (phaseType == "SUBI") {
        // Ionic sublattice
        // computeExcessGibbsSUBI(ctx, phaseIndex);
    }
}

// Wrapper function for compatibility with Subminimization
// This is the function signature expected by the Fortran-style code
void compExcessGibbsEnergy(ThermoContext& ctx, int phaseIndex) {
    // NOTE: Subminimization passes 1-based index, but computeExcessGibbs expects 0-based
    // So we subtract 1 here to convert
    computeExcessGibbs(ctx, phaseIndex - 1);
}

} // namespace Thermochimica
