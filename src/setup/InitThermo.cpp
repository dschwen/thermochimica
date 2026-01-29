#include "thermochimica/Thermochimica.hpp"

namespace Thermochimica {

// Additional initialization routines

void initTolerances(ThermoContext& ctx) {
    ctx.thermo->tolerances.initDefaults();
}

} // namespace Thermochimica
