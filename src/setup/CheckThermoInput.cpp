#include "thermochimica/Thermochimica.hpp"
#include "thermochimica/util/ErrorCodes.hpp"

namespace Thermochimica {

// Input validation functions - called by checkSystem()

bool validateTemperature(const ThermoContext& ctx) {
    double T = ctx.io->dTemperature;
    return T > 0.0 && T <= 10000.0;
}

bool validatePressure(const ThermoContext& ctx) {
    double P = ctx.io->dPressure;
    return P > 0.0 && P <= 1e10;
}

bool validateComposition(const ThermoContext& ctx) {
    for (int i = 0; i < Constants::kMaxIsotopes; ++i) {
        if (ctx.io->dElementMass[i] > 0.0) {
            return true;
        }
    }
    return false;
}

} // namespace Thermochimica
