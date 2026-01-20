#include "thermochimica/ThermoContext.hpp"

namespace Thermochimica {

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
        // computeExcessGibbsQKTO(ctx, phaseIndex);
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

} // namespace Thermochimica
