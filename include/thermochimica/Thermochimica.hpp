#pragma once

#include "ThermoContext.hpp"
#include <string>

namespace Thermochimica {

// ============================================================================
// IMPORTANT: Use ThermoClass for all thermodynamic calculations
// ============================================================================
// The free-function API has been removed in favor of the modern object-oriented
// ThermoClass API. Please use ThermoClass instead.
//
// Migration example:
//   Old (removed):
//     ThermoContext ctx;
//     setThermoFilename(ctx, "CO.dat");
//     thermochimica(ctx);
//
//   New (recommended):
//     ThermoClass thermo;
//     thermo.loadDatabase("CO.dat");
//     thermo.calculate();
//
// See ThermoClass.hpp for the complete API.
// ============================================================================

// ============================================================================
// Utility Functions
// ============================================================================

/// Get atomic number from element symbol
int getAtomicNumber(const std::string& symbol);

/// Get element symbol from atomic number
std::string getElementSymbol(int atomicNumber);

/// Get error message for error code
const char* getErrorMessage(int errorCode);

} // namespace Thermochimica
