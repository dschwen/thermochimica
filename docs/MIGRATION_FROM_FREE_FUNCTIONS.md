# Migration Guide: Free-Function API → ThermoClass

[← Back to README](../README.md) | [ThermoClass API Reference](class_based_api.md)

---

## Overview

The free-function API (where you pass `ThermoContext&` to functions) has been **removed** in favor of the modern object-oriented `ThermoClass` API. This guide helps you migrate existing code to the new API.

## Why the Change?

The class-based API provides significant advantages:

- **Modern C++ Design**: Natural object-oriented interface using RAII and encapsulation
- **Type Safety**: Member functions prevent API misuse and improve compile-time checking
- **Better IDE Support**: Autocomplete and discoverability work much better with member functions
- **Strategy Pattern**: Easy injection of custom solvers, models, and line search algorithms
- **Move Semantics**: Efficient resource management with move constructors
- **Cleaner Code**: No need to pass context to every function call
- **Identical Performance**: Zero overhead compared to free functions (all tests produce identical results)

## Quick Migration Example

### Before (Free-Function API)

```cpp
#include <thermochimica/Thermochimica.hpp>

int main() {
    Thermochimica::ThermoContext ctx;

    Thermochimica::setStandardUnits(ctx);
    Thermochimica::setThermoFilename(ctx, "CO.dat");
    Thermochimica::parseCSDataFile(ctx);

    Thermochimica::setTemperaturePressure(ctx, 1000.0, 1.0);
    Thermochimica::setElementMass(ctx, 6, 1.0);
    Thermochimica::setElementMass(ctx, 8, 2.0);

    Thermochimica::thermochimica(ctx);

    if (ctx.isSuccess()) {
        double G = Thermochimica::getGibbsEnergy(ctx);
        auto [moles, info] = Thermochimica::getMolesPhase(ctx, "gas_ideal");
        Thermochimica::printResults(ctx);
    }

    return 0;
}
```

### After (ThermoClass API)

```cpp
#include <thermochimica/ThermoClass.hpp>

int main() {
    Thermochimica::ThermoClass thermo;

    thermo.setStandardUnits();
    thermo.loadDatabase("CO.dat");  // Combined setThermoFilename + parseCSDataFile

    thermo.setTemperaturePressure(1000.0, 1.0);
    thermo.setElementMass(6, 1.0);
    thermo.setElementMass(8, 2.0);

    thermo.calculate();  // Was: thermochimica(ctx)

    if (thermo.isSuccess()) {
        double G = thermo.getGibbsEnergy();
        auto [moles, info] = thermo.getMolesPhase("gas_ideal");
        thermo.printResults();
    }

    return 0;
}
```

## Migration Steps

1. **Change include**: `#include <thermochimica/ThermoClass.hpp>`
2. **Replace ThermoContext**: `ThermoContext ctx` → `ThermoClass thermo`
3. **Remove `ctx,` argument**: All function calls become member function calls
4. **Update function names**:
   - `setThermoFilename() + parseCSDataFile()` → `loadDatabase()`
   - `thermochimica()` → `calculate()`
5. **Update method calls**: `Thermochimica::func(ctx, args)` → `thermo.func(args)`

## Complete API Mapping

### Database Loading

| Old (Free-Function) | New (ThermoClass) |
|---------------------|-------------------|
| `setThermoFilename(ctx, "file.dat")` | `setThermoFilename("file.dat")` |
| `parseCSDataFile(ctx)` | `parseCSDataFile()` |
| `setThermoFilename(ctx, f); parseCSDataFile(ctx);` | `loadDatabase("file.dat")` |

### Input Configuration

| Old (Free-Function) | New (ThermoClass) |
|---------------------|-------------------|
| `setTemperature(ctx, T)` | `setTemperature(T)` |
| `setPressure(ctx, P)` | `setPressure(P)` |
| `setTemperaturePressure(ctx, T, P)` | `setTemperaturePressure(T, P)` |
| `setElementMass(ctx, Z, mass)` | `setElementMass(Z, mass)` |
| `setElementMass(ctx, "C", mass)` | `setElementMass("C", mass)` |
| `setStandardUnits(ctx)` | `setStandardUnits()` |
| `setUnitsSI(ctx)` | `setUnitsSI()` |
| `setUnitTemperature(ctx, unit)` | `setUnitTemperature(unit)` |
| `setUnitPressure(ctx, unit)` | `setUnitPressure(unit)` |
| `setUnitMass(ctx, unit)` | `setUnitMass(unit)` |

### Calculation

| Old (Free-Function) | New (ThermoClass) |
|---------------------|-------------------|
| `thermochimica(ctx)` | `calculate()` |
| `init(ctx)` | `initialize()` |
| `checkSystem(ctx)` | `checkSystem()` |
| `compThermoData(ctx)` | `computeThermoData()` |
| `setup(ctx)` | `setup()` |
| `solve(ctx)` | `solve()` |

### Phase Constraints

| Old (Free-Function) | New (ThermoClass) |
|---------------------|-------------------|
| `setSolnPhaseConstraint(ctx, name, frac)` | `setSolnPhaseConstraint(name, frac)` |
| `setCondPhaseConstraint(ctx, name, frac)` | `setCondPhaseConstraint(name, frac)` |
| `removePhaseConstraint(ctx, name)` | `removePhaseConstraint(name)` |
| `clearPhaseConstraints(ctx)` | `clearPhaseConstraints()` |
| `getPhaseElementFraction(ctx, name)` | `getPhaseElementFraction(name)` |
| `arePhaseConstraintsSatisfied(ctx)` | `arePhaseConstraintsSatisfied()` |

### Output Retrieval

| Old (Free-Function) | New (ThermoClass) |
|---------------------|-------------------|
| `getGibbsEnergy(ctx)` | `getGibbsEnergy()` |
| `getHeatCapacity(ctx)` | `getHeatCapacity()` |
| `getEntropy(ctx)` | `getEntropy()` |
| `getEnthalpy(ctx)` | `getEnthalpy()` |
| `getMolesPhase(ctx, name)` | `getMolesPhase(name)` |
| `getOutputChemPot(ctx, elem)` | `getOutputChemPot(elem)` |
| `getOutputSolnSpecies(ctx, ph, sp)` | `getOutputSolnSpecies(ph, sp)` |
| `getOutputMolSpecies(ctx, ph, sp)` | `getOutputMolSpecies(ph, sp)` |
| `printResults(ctx)` | `printResults()` |

### Status and Error Handling

| Old (Free-Function) | New (ThermoClass) |
|---------------------|-------------------|
| `ctx.infoThermo()` | `getInfoCode()` |
| `ctx.isSuccess()` | `isSuccess()` |
| `getErrorMessage(ctx.infoThermo())` | `getErrorMessage()` |

### Reset Functions

| Old (Free-Function) | New (ThermoClass) |
|---------------------|-------------------|
| `resetThermo(ctx)` | `reset()` |
| `resetThermoAll(ctx)` | `resetAll()` |

### Database Queries

| Old (Free-Function) | New (ThermoClass) |
|---------------------|-------------------|
| `getNumberElementsDatabase(ctx)` | `getNumberElementsDatabase()` |
| `getElementAtIndex(ctx, i)` | `getElementAtIndex(i)` |
| `getNumberPhasesSystem(ctx)` | `getNumberPhasesSystem()` |
| `getPhaseNameAtIndex(ctx, i)` | `getPhaseNameAtIndex(i)` |
| `getNumberSpeciesSystem(ctx)` | `getNumberSpeciesSystem()` |
| `getSpeciesAtIndex(ctx, i)` | `getSpeciesAtIndex(i)` |

## Common Migration Patterns

### Pattern 1: Basic Calculation

**Before:**
```cpp
Thermochimica::ThermoContext ctx;
Thermochimica::setStandardUnits(ctx);
Thermochimica::setThermoFilename(ctx, "CO.dat");
Thermochimica::parseCSDataFile(ctx);
Thermochimica::setTemperaturePressure(ctx, 1000.0, 1.0);
Thermochimica::setElementMass(ctx, 6, 1.0);
Thermochimica::thermochimica(ctx);
double G = Thermochimica::getGibbsEnergy(ctx);
```

**After:**
```cpp
Thermochimica::ThermoClass thermo;
thermo.setStandardUnits();
thermo.loadDatabase("CO.dat");
thermo.setTemperaturePressure(1000.0, 1.0);
thermo.setElementMass(6, 1.0);
thermo.calculate();
double G = thermo.getGibbsEnergy();
```

### Pattern 2: Multi-Point Calculations

**Before:**
```cpp
Thermochimica::ThermoContext ctx;
Thermochimica::setThermoFilename(ctx, "CO.dat");
Thermochimica::parseCSDataFile(ctx);
Thermochimica::setStandardUnits(ctx);

for (double T = 500; T <= 2000; T += 100) {
    Thermochimica::resetThermo(ctx);
    Thermochimica::setTemperaturePressure(ctx, T, 1.0);
    Thermochimica::setElementMass(ctx, 6, 1.0);
    Thermochimica::setElementMass(ctx, 8, 2.0);
    Thermochimica::thermochimica(ctx);
    if (ctx.isSuccess()) {
        double G = Thermochimica::getGibbsEnergy(ctx);
        std::cout << "T=" << T << ", G=" << G << "\n";
    }
}
```

**After:**
```cpp
Thermochimica::ThermoClass thermo;
thermo.loadDatabase("CO.dat");
thermo.setStandardUnits();

for (double T = 500; T <= 2000; T += 100) {
    thermo.reset();
    thermo.setTemperaturePressure(T, 1.0);
    thermo.setElementMass(6, 1.0);
    thermo.setElementMass(8, 2.0);
    thermo.calculate();
    if (thermo.isSuccess()) {
        double G = thermo.getGibbsEnergy();
        std::cout << "T=" << T << ", G=" << G << "\n";
    }
}
```

### Pattern 3: Thread-Safe Parallel Calculations

**Before:**
```cpp
void calculate(double T) {
    Thermochimica::ThermoContext ctx;
    Thermochimica::setThermoFilename(ctx, "CO.dat");
    Thermochimica::parseCSDataFile(ctx);
    Thermochimica::setStandardUnits(ctx);
    Thermochimica::setTemperaturePressure(ctx, T, 1.0);
    Thermochimica::setElementMass(ctx, 6, 1.0);
    Thermochimica::thermochimica(ctx);
}

std::vector<std::thread> threads;
for (double T = 500; T <= 2000; T += 100) {
    threads.emplace_back(calculate, T);
}
```

**After:**
```cpp
void calculate(double T) {
    Thermochimica::ThermoClass thermo;
    thermo.loadDatabase("CO.dat");
    thermo.setStandardUnits();
    thermo.setTemperaturePressure(T, 1.0);
    thermo.setElementMass(6, 1.0);
    thermo.calculate();
}

std::vector<std::thread> threads;
for (double T = 500; T <= 2000; T += 100) {
    threads.emplace_back(calculate, T);
}
```

### Pattern 4: Phase Constraints

**Before:**
```cpp
Thermochimica::ThermoContext ctx;
Thermochimica::setThermoFilename(ctx, "NobleMetals-Kaye.dat");
Thermochimica::parseCSDataFile(ctx);
Thermochimica::setStandardUnits(ctx);
Thermochimica::setTemperaturePressure(ctx, 1500.0, 1.0);
Thermochimica::setElementMass(ctx, 42, 0.5);
Thermochimica::setElementMass(ctx, 44, 0.5);
Thermochimica::setSolnPhaseConstraint(ctx, "FCCN", 0.4);
Thermochimica::setSolnPhaseConstraint(ctx, "BCCN", 0.6);
Thermochimica::thermochimica(ctx);
auto [frac, info] = Thermochimica::getPhaseElementFraction(ctx, "FCCN");
```

**After:**
```cpp
Thermochimica::ThermoClass thermo;
thermo.loadDatabase("NobleMetals-Kaye.dat");
thermo.setStandardUnits();
thermo.setTemperaturePressure(1500.0, 1.0);
thermo.setElementMass(42, 0.5);
thermo.setElementMass(44, 0.5);
thermo.setSolnPhaseConstraint("FCCN", 0.4);
thermo.setSolnPhaseConstraint("BCCN", 0.6);
thermo.calculate();
auto [frac, info] = thermo.getPhaseElementFraction("FCCN");
```

### Pattern 5: Accessing Internal State

**Before:**
```cpp
Thermochimica::ThermoContext ctx;
// ... setup and calculate ...
auto& thermo = *ctx.thermo;
for (int i = 0; i < thermo.nSpecies; ++i) {
    double x = thermo.dMolFraction(i);
}
```

**After:**
```cpp
Thermochimica::ThermoClass thermo;
// ... setup and calculate ...
ThermoContext& ctx = thermo.getContext();
auto& state = *ctx.thermo;
for (int i = 0; i < state.nSpecies; ++i) {
    double x = state.dMolFraction(i);
}
```

### Pattern 6: Error Handling

**Before:**
```cpp
Thermochimica::ThermoContext ctx;
Thermochimica::setThermoFilename(ctx, "file.dat");
Thermochimica::parseCSDataFile(ctx);
if (ctx.infoThermo() != 0) {
    std::cerr << "Parse failed: "
              << Thermochimica::getErrorMessage(ctx.infoThermo()) << "\n";
}

Thermochimica::thermochimica(ctx);
if (!ctx.isSuccess()) {
    std::cerr << "Calculation failed\n";
}
```

**After:**
```cpp
Thermochimica::ThermoClass thermo;
int result = thermo.loadDatabase("file.dat");
if (result != 0) {
    std::cerr << "Parse failed: " << thermo.getErrorMessage() << "\n";
}

thermo.calculate();
if (!thermo.isSuccess()) {
    std::cerr << "Calculation failed: " << thermo.getErrorMessage() << "\n";
}
```

## Breaking Changes

### 1. Include Path Changed

**Old:** `#include <thermochimica/Thermochimica.hpp>`
**New:** `#include <thermochimica/ThermoClass.hpp>`

### 2. Main Class Name Changed

**Old:** `Thermochimica::ThermoContext ctx;`
**New:** `Thermochimica::ThermoClass thermo;`

### 3. Function Call Syntax Changed

**Old:** `Thermochimica::function(ctx, args...)`
**New:** `thermo.function(args...)`

### 4. Main Solver Function Renamed

**Old:** `Thermochimica::thermochimica(ctx);`
**New:** `thermo.calculate();`

### 5. Combined Database Loading

**Old:**
```cpp
Thermochimica::setThermoFilename(ctx, "file.dat");
Thermochimica::parseCSDataFile(ctx);
```

**New:**
```cpp
thermo.loadDatabase("file.dat");
// Or separately if needed:
thermo.setThermoFilename("file.dat");
thermo.parseCSDataFile();
```

### 6. Status Checking Methods

**Old:** `ctx.infoThermo()`, `ctx.isSuccess()`
**New:** `thermo.getInfoCode()`, `thermo.isSuccess()`

### 7. Error Messages

**Old:** `Thermochimica::getErrorMessage(ctx.infoThermo())`
**New:** `thermo.getErrorMessage()`

### 8. Namespace Changes

All API functions remain in `Thermochimica::` namespace, but are now member functions of `ThermoClass`.

## Advanced Features

### Custom Solvers (New in ThermoClass)

ThermoClass supports dependency injection for custom solvers:

```cpp
#include <thermochimica/interfaces/ISolver.hpp>

class MyCustomSolver : public ISolver {
    // Implement custom solver logic
};

Thermochimica::ThermoClass thermo;
thermo.setSolver(std::make_unique<MyCustomSolver>());
thermo.calculate();  // Uses custom solver
```

This was not possible with the free-function API.

### Move Semantics (New in ThermoClass)

```cpp
Thermochimica::ThermoClass createThermo() {
    Thermochimica::ThermoClass thermo;
    thermo.loadDatabase("CO.dat");
    thermo.setStandardUnits();
    return thermo;  // Efficient move, not copy
}

Thermochimica::ThermoClass thermo = createThermo();
thermo.calculate();
```

### RAII Benefits (New in ThermoClass)

Resources are automatically managed:

```cpp
{
    Thermochimica::ThermoClass thermo;
    thermo.loadDatabase("CO.dat");
    thermo.calculate();
    // Automatic cleanup on scope exit
}
```

## Interoperability with Legacy Code

If you need to gradually migrate or interoperate with legacy code using the old API:

```cpp
#include <thermochimica/ThermoClass.hpp>
#include <thermochimica/Thermochimica.hpp>

// New code using ThermoClass
Thermochimica::ThermoClass thermo;
thermo.loadDatabase("CO.dat");
thermo.setStandardUnits();
thermo.setTemperaturePressure(1000.0, 1.0);

// Access underlying context for legacy functions
ThermoContext& ctx = thermo.getContext();

// Can use legacy free functions if needed
Thermochimica::setElementMass(ctx, 6, 1.0);  // Legacy call

// Back to new API
thermo.calculate();
```

**Note:** This is a temporary migration strategy. The free-function API has been removed, so you should fully migrate to ThermoClass.

## Troubleshooting

### Issue: "ThermoContext not found"

**Cause:** Using old include
**Solution:** Change to `#include <thermochimica/ThermoClass.hpp>`

### Issue: "No matching function for call to 'thermochimica'"

**Cause:** Old function name
**Solution:** Change `thermochimica(ctx)` to `thermo.calculate()`

### Issue: "Too many arguments to function"

**Cause:** Still passing context as first argument
**Solution:** Remove `ctx,` from function calls

### Issue: "ctx not declared in this scope"

**Cause:** Code still references old context variable
**Solution:** Replace `ctx` with `thermo` and use member functions

### Issue: Accessing internal state

**Cause:** Need access to ThermoContext for advanced use
**Solution:** Use `thermo.getContext()` to access underlying context

## Performance Considerations

The ThermoClass API has **zero performance overhead** compared to the free-function API:

- Virtual function calls: < 1% overhead (only one virtual call per strategy component)
- Memory overhead: Minimal (stores strategy pointers)
- All 90 tests produce numerically identical results
- Same convergence behavior and iteration counts

## Benefits Summary

| Aspect | Free-Function API | ThermoClass API |
|--------|-------------------|-----------------|
| Syntax | `func(ctx, args)` | `thermo.func(args)` |
| Encapsulation | Manual context passing | Automatic |
| Type Safety | Moderate | High |
| IDE Support | Limited | Excellent |
| Extensibility | Difficult | Easy (strategy pattern) |
| Modern C++ | Partial | Full (RAII, move) |
| Performance | Baseline | Identical |

## Getting Help

- **API Documentation**: See [class_based_api.md](class_based_api.md)
- **Examples**: See [examples.md](examples.md) for updated examples
- **Issues**: File issues on GitHub with migration questions

---

[← Back to README](../README.md) | [ThermoClass API Reference](class_based_api.md)
