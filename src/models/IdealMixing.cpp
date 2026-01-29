#include "thermochimica/ThermoContext.hpp"
#include <cmath>

namespace Thermochimica {

// Ideal mixing model (IDMX)
// G^xs = 0, only configurational entropy contribution

void computeExcessGibbsIDMX(ThermoContext& ctx, int phaseIndex) {
    auto& thermo = *ctx.thermo;

    // For ideal mixing, excess Gibbs energy is zero
    thermo.dGibbsSolnPhase(phaseIndex) = 0.0;

    // Compute ideal mixing contribution to chemical potential
    // mu_i = mu_i^0 + RT*ln(x_i)

    int iFirst = (phaseIndex > 0) ? thermo.nSpeciesPhase(phaseIndex) : 0;
    int iLast = thermo.nSpeciesPhase(phaseIndex + 1);

    for (int i = iFirst; i < iLast; ++i) {
        double x = thermo.dMolFraction(i);
        if (x > 1e-100) {
            // Add RT*ln(x) contribution (already normalized by RT)
            thermo.dChemicalPotential(i) += std::log(x);
        }
    }
}

} // namespace Thermochimica
