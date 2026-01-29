# Thermochimica C++ API Migration Guide

This guide helps users migrate from the Fortran-based C++ API to the new pure C++ implementation.

## Overview

The main difference between the old and new APIs is **state management**:

| Aspect | Old Fortran API | New C++ API |
|--------|-----------------|-------------|
| State | Global (hidden in Fortran modules) | Explicit `ThermoContext` object |
| Thread safety | Not thread-safe | Thread-safe (each thread uses own context) |
| Multiple systems | Not supported | Supported (multiple contexts) |
| Memory management | Automatic (Fortran) | RAII via smart pointers |

## Quick Start

### Old API (Fortran-backed)
```cpp
#include "Thermochimica-cxx.h"

int main() {
    Thermochimica::setThermoFilename("database.dat");
    Thermochimica::parseThermoFile();

    Thermochimica::setStandardUnits();
    Thermochimica::setTemperaturePressure(1000.0, 1.0);
    Thermochimica::setElementMass(26, 0.5);  // Fe
    Thermochimica::setElementMass(8, 0.5);   // O

    Thermochimica::thermochimica();

    if (Thermochimica::checkInfoThermo() == 0) {
        auto [gibbs, info] = Thermochimica::getSolnPhaseMol("FCC_A1");
        // ...
    }

    Thermochimica::resetThermo();
    return 0;
}
```

### New API
```cpp
#include <thermochimica/Thermochimica.hpp>

int main() {
    Thermochimica::ThermoContext ctx;

    Thermochimica::setThermoFilename(ctx, "database.dat");
    Thermochimica::parseCSDataFile(ctx);

    Thermochimica::setStandardUnits(ctx);
    Thermochimica::setTemperaturePressure(ctx, 1000.0, 1.0);
    Thermochimica::setElementMass(ctx, 26, 0.5);  // Fe
    Thermochimica::setElementMass(ctx, 8, 0.5);   // O

    Thermochimica::thermochimica(ctx);

    if (ctx.infoThermo() == 0) {
        auto [moles, info] = Thermochimica::getSolnPhaseMol(ctx, "FCC_A1");
        // ...
    }

    Thermochimica::resetThermo(ctx);
    return 0;
}
```

## Function-by-Function Migration

### Setup Functions

| Old API | New API |
|---------|---------|
| `setThermoFilename(filename)` | `setThermoFilename(ctx, filename)` |
| `parseThermoFile()` | `parseCSDataFile(ctx)` |
| `setStandardUnits()` | `setStandardUnits(ctx)` |
| `setModelicaUnits()` | `setUnits(ctx, "K", "Pa", "kg/m3")` |
| `setUnitTemperature(unit)` | `setUnitTemperature(ctx, unit)` |
| `setUnitPressure(unit)` | `setUnitPressure(ctx, unit)` |
| `setUnitMass(unit)` | `setUnitMass(ctx, unit)` |
| `setTemperaturePressure(T, P)` | `setTemperaturePressure(ctx, T, P)` |
| `setElementMass(Z, mass)` | `setElementMass(ctx, Z, mass)` |
| `presetElementMass(Z, mass)` | `setElementMass(ctx, Z, mass)` |

### Solver Functions

| Old API | New API |
|---------|---------|
| `thermochimica()` | `thermochimica(ctx)` |
| `setup()` | `setup(ctx)` |
| `solve()` | `solve(ctx)` |
| `init()` | `init(ctx)` |
| `checkSystem()` | `checkSystem(ctx)` |
| `compThermoData()` | `compThermoData(ctx)` |

### Status and Reset Functions

| Old API | New API |
|---------|---------|
| `checkInfoThermo()` | `ctx.infoThermo()` |
| `resetInfoThermo()` | `ctx.io->INFOThermo = 0` |
| `resetThermo()` | `resetThermo(ctx)` |
| `resetThermoAll()` | `resetThermoAll(ctx)` |

### Output Functions

| Old API | New API |
|---------|---------|
| `printResults()` | `printResults(ctx)` |
| `setPrintResultsMode(mode)` | `setPrintResultsMode(ctx, mode)` |

### Data Retrieval Functions

| Old API | New API |
|---------|---------|
| `getOutputChemPot(elem)` | `getOutputChemPot(ctx, elem)` |
| `getOutputSolnSpecies(phase, species)` | `getOutputSolnSpecies(ctx, phase, species)` |
| `getOutputMolSpecies(species)` | `getOutputMolSpecies(ctx, phase, species)` |
| `getSolnPhaseMol(phase)` | `getSolnPhaseMol(ctx, phase)` |
| `getPureConPhaseMol(phase)` | `getPureConPhaseMol(ctx, phase)` |
| `getElementMolesInPhase(elem, phase)` | `getElementMolesInPhase(ctx, elem, phase)` |
| `getPhaseIndex(phase)` | `getPhaseIndex(ctx, phase)` |
| `getOutputSiteFraction(phase, sub, const)` | `getOutputSiteFraction(ctx, phase, sub, constName)` |
| `isPhaseGas(index)` | `isPhaseGas(ctx, phaseName)` |
| `isPhaseMQM(index)` | `isPhaseMQM(ctx, phaseName)` |

### Database Query Functions

| Old API | New API |
|---------|---------|
| `getNumberElementsDatabase()` | `getNumberElementsDatabase(ctx)` |
| `getElementsDatabase()` | Loop with `getElementAtIndex(ctx, i)` |
| `getElementAtIndex(index)` | `getElementAtIndex(ctx, index)` |
| `getNumberPhasesSystem()` | `getNumberPhasesSystem(ctx)` |
| `getPhaseNamesSystem()` | Loop with `getPhaseNameAtIndex(ctx, i)` |
| `getPhaseNameAtIndex(index)` | `getPhaseNameAtIndex(ctx, index)` |
| `getNumberSpeciesSystem()` | `getNumberSpeciesSystem(ctx)` |
| `getSpeciesInPhase(index)` | Loop with `getSpeciesAtIndex(ctx, i)` |

### Thermodynamic Properties

| Old API | New API |
|---------|---------|
| `setHeatCapacityEnthalpyEntropyRequested(b)` | `setHeatCapacityEntropyEnthalpy(ctx, b)` |
| `getHeatCapacityEnthalpyEntropy()` | Individual: `getHeatCapacity(ctx)`, `getEnthalpy(ctx)`, `getEntropy(ctx)` |
| N/A | `getGibbsEnergy(ctx)` |

### Reinitialization Functions

| Old API | New API |
|---------|---------|
| `saveReinitData()` | `saveReinitData(ctx)` |
| `getReinitDataSizes()` | N/A (handled internally) |
| `getReinitData()` | `getReinitData(ctx)` |
| `setReinitData(data)` | `setReinitData(ctx, assemblage, chemPot, molesPhase, elemPot, molFrac)` |
| `setReinitRequested(b)` | `setReinitRequested(ctx, b)` |
| `resetReinit()` | Part of `resetThermo(ctx)` |

### Advanced Functions

| Old API | New API |
|---------|---------|
| `setFuzzyStoich(b)` | `setFuzzyStoich(ctx, b)` |
| `setFuzzyMagnitude(mag)` | `setFuzzyMagnitude(ctx, mag)` |
| `setGibbsMinCheck(b)` | N/A (always performed) |
| `getMqmqaMolesPairs(phase)` | To be implemented |
| `getMqmqaPairMolFraction(phase, pair)` | To be implemented |
| `getMqmqaNumberPairsQuads(phase)` | To be implemented |
| `getMqmqaConstituentFraction(phase, sub, const)` | To be implemented |

## Key Differences

### 1. Context Object

Every function now takes a `ThermoContext&` as the first parameter. This object contains all state:

```cpp
Thermochimica::ThermoContext ctx;

// ctx contains:
// - ctx.thermo: Thermodynamic state (species, phases, results)
// - ctx.io: Input/output parameters (T, P, composition, output flags)
// - ctx.gem: GEM solver state
// - ctx.parser: Parser state (database coefficients)
// - ctx.reinit: Reinitialization state
```

### 2. Thread Safety

The new API is thread-safe because each thread can have its own context:

```cpp
void calculateEquilibrium(int threadId, double temp) {
    Thermochimica::ThermoContext ctx;  // Thread-local context

    Thermochimica::setThermoFilename(ctx, "database.dat");
    Thermochimica::parseCSDataFile(ctx);
    Thermochimica::setStandardUnits(ctx);
    Thermochimica::setTemperaturePressure(ctx, temp, 1.0);
    // ... set composition ...
    Thermochimica::thermochimica(ctx);

    // Results are in ctx, no global state conflicts
}
```

### 3. Multiple Simultaneous Systems

You can now work with multiple thermodynamic systems simultaneously:

```cpp
Thermochimica::ThermoContext ctx1, ctx2;

// Load different databases
Thermochimica::setThermoFilename(ctx1, "steel.dat");
Thermochimica::setThermoFilename(ctx2, "slag.dat");

Thermochimica::parseCSDataFile(ctx1);
Thermochimica::parseCSDataFile(ctx2);

// Calculate independently
Thermochimica::thermochimica(ctx1);
Thermochimica::thermochimica(ctx2);
```

### 4. Error Handling

```cpp
// Old API
Thermochimica::thermochimica();
int info = Thermochimica::checkInfoThermo();
if (info != 0) {
    // Error occurred
}

// New API
Thermochimica::thermochimica(ctx);
if (ctx.infoThermo() != 0) {
    const char* msg = Thermochimica::getErrorMessage(ctx.infoThermo());
    std::cerr << "Error: " << msg << std::endl;
}
```

### 5. Return Types

Some functions have changed return types:

```cpp
// Old: Returns pair<double, int>
auto [moles, info] = Thermochimica::getSolnPhaseMol("FCC_A1");

// New: Same return type, but requires context
auto [moles, info] = Thermochimica::getSolnPhaseMol(ctx, "FCC_A1");
```

### 6. Phase/Species Lookup

The new API prefers string names over indices for clarity:

```cpp
// Old API used index for some functions
bool isGas = Thermochimica::isPhaseGas(0);  // Phase index 0

// New API uses phase name
bool isGas = Thermochimica::isPhaseGas(ctx, "gas_ideal");
```

## Compatibility Wrapper (Optional)

For minimal code changes during migration, you can create a compatibility wrapper:

```cpp
// ThermochimicaCompat.hpp - Compatibility layer
#pragma once
#include <thermochimica/Thermochimica.hpp>

namespace Thermochimica {
namespace Compat {

// Global context for compatibility
inline ThermoContext& globalContext() {
    static ThermoContext ctx;
    return ctx;
}

// Wrapper functions matching old API signatures
inline void setThermoFilename(const std::string& f) {
    Thermochimica::setThermoFilename(globalContext(), f);
}

inline void parseThermoFile() {
    Thermochimica::parseCSDataFile(globalContext());
}

inline void setStandardUnits() {
    Thermochimica::setStandardUnits(globalContext());
}

inline void setTemperaturePressure(double T, double P) {
    Thermochimica::setTemperaturePressure(globalContext(), T, P);
}

inline void setElementMass(int Z, double mass) {
    Thermochimica::setElementMass(globalContext(), Z, mass);
}

inline void thermochimica() {
    Thermochimica::thermochimica(globalContext());
}

inline int checkInfoThermo() {
    return globalContext().infoThermo();
}

inline std::pair<double, int> getSolnPhaseMol(const std::string& phase) {
    return Thermochimica::getSolnPhaseMol(globalContext(), phase);
}

inline void resetThermo() {
    Thermochimica::resetThermo(globalContext());
}

// Add more wrappers as needed...

}} // namespace Thermochimica::Compat
```

Usage with compatibility wrapper:

```cpp
#include "ThermochimicaCompat.hpp"
using namespace Thermochimica::Compat;

int main() {
    // Code looks almost identical to old API
    setThermoFilename("database.dat");
    parseThermoFile();
    setStandardUnits();
    setTemperaturePressure(1000.0, 1.0);
    setElementMass(26, 0.5);
    thermochimica();

    if (checkInfoThermo() == 0) {
        auto [moles, info] = getSolnPhaseMol("FCC_A1");
    }
    return 0;
}
```

## Recommended Migration Strategy

1. **Phase 1**: Add `#include <thermochimica/Thermochimica.hpp>` and create a single `ThermoContext ctx` at the start of your calculation routine.

2. **Phase 2**: Add `ctx` as the first argument to all Thermochimica function calls.

3. **Phase 3**: Replace `checkInfoThermo()` with `ctx.infoThermo()`.

4. **Phase 4**: Update any index-based phase lookups to use phase names.

5. **Phase 5**: Test thoroughly, especially reinitialization and multi-step calculations.

## Benefits of Migration

- **Thread safety**: Run multiple calculations in parallel
- **Multiple systems**: Work with different databases simultaneously
- **Better debugging**: All state is visible in the context object
- **No global state**: Easier to test and reason about
- **Modern C++**: RAII, smart pointers, standard containers

## Building

The new C++ library uses CMake:

```bash
cd cpp
mkdir build && cd build
cmake ..
make -j4

# Link against libthermochimica.a
# Include path: cpp/include
```

## Questions?

If you encounter issues during migration, please open an issue on GitHub.
