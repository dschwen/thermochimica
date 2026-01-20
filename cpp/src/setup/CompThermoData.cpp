#include "thermochimica/Thermochimica.hpp"
#include <cmath>

namespace Thermochimica {

// Compute thermodynamic data - implementation details

void computeGibbsEnergy(ThermoContext& ctx, int speciesIndex, double T) {
    // Compute Gibbs energy from polynomial coefficients
    // G = A1 + A2*T + A3*T*ln(T) + A4*T^2 + A5*T^3 + A6/T + ...

    // This would use the Gibbs coefficients parsed from the data file
    // Placeholder implementation
}

void computeChemicalPotentials(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;
    double T = ctx.io->dTemperature;
    double R = Constants::kIdealGasConstant;

    for (int i = 0; i < thermo.nSpecies; ++i) {
        // mu_i / RT = G_i / RT
        thermo.dChemicalPotential(i) = thermo.dStdGibbsEnergy(i) / (R * T);
    }
}

} // namespace Thermochimica
