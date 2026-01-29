# Error Handling

[← Back to README](../README.md) | [Examples](examples.md) | [Architecture →](architecture.md)

---

How to handle errors and troubleshoot common problems.

## Checking for Errors

After any operation, check the error code:

```cpp
Thermochimica::thermochimica(ctx);

if (ctx.infoThermo() != 0) {
    std::cerr << "Error " << ctx.infoThermo() << ": "
              << Thermochimica::getErrorMessage(ctx.infoThermo())
              << std::endl;
}

// Or use the convenience method
if (!ctx.isSuccess()) {
    // Handle error
}
```

---

## Error Code Reference

### Success (0)

| Code | Name | Description |
|------|------|-------------|
| 0 | `kSuccess` | Operation completed successfully |

### Input Validation Errors (1-9)

| Code | Name | Description |
|------|------|-------------|
| 1 | `kTemperatureOutOfRange` | Temperature is negative, zero, or invalid |
| 2 | `kPressureOutOfRange` | Pressure is negative, zero, or invalid |
| 3 | `kCompositionOutOfRange` | Negative mass or zero total composition |
| 4 | `kUnrecognizedTemperatureUnit` | Invalid temperature unit string |
| 5 | `kUnrecognizedPressureUnit` | Invalid pressure unit string |
| 6 | `kNoDataFileSpecified` | No database filename set before calculation |
| 7 | `kUnrecognizedMassUnit` | Invalid mass unit string |
| 8 | `kInvalidCompoundStoichiometry` | Invalid compound stoichiometry specified |

### Parser Errors (10-19)

| Code | Name | Description |
|------|------|-------------|
| 10 | `kDataFileNotFound` | Database file does not exist |
| 11 | `kParserError` | Error parsing database file format |
| 12 | `kUnsupportedPhaseType` | Database contains unsupported phase model |
| 13 | `kInvalidDataFormat` | Invalid data format in database file |

### System Check Errors (20-29)

| Code | Name | Description |
|------|------|-------------|
| 20 | `kNoElementsInSystem` | No elements defined in system |
| 21 | `kNoSpeciesInSystem` | No species found for given elements |
| 22 | `kSystemCheckFailed` | General system validation failure |
| 23 | `kInsufficientElements` | Not enough elements for calculation |

### Solver Errors (30-49)

| Code | Name | Description |
|------|------|-------------|
| 30 | `kGEMSolverDidNotConverge` | GEM solver exceeded max iterations |
| 31 | `kLevelingSolverFailed` | Initial leveling solver failed |
| 32 | `kSingularMatrix` | Singular matrix in linear solve |
| 33 | `kNoStablePhases` | No thermodynamically stable phases |
| 34 | `kMaxIterationsReached` | Maximum iterations reached |
| 35 | `kNumericalInstability` | Numerical instability detected |
| 36 | `kPhaseRuleViolation` | Gibbs phase rule violated |

### Post-Processing Errors (50-59)

| Code | Name | Description |
|------|------|-------------|
| 50 | `kPostProcessError` | Error in post-processing calculations |
| 51 | `kOutputError` | Error writing output |
| 52 | `kJSONWriteError` | Error writing JSON file |

### Reinitialization Errors (60-69)

| Code | Name | Description |
|------|------|-------------|
| 60 | `kReinitDataNotAvailable` | No saved reinitialization data |
| 61 | `kReinitDataInvalid` | Reinitialization data incompatible |

---

## Common Problems and Solutions

### "Data file not found" (Error 10)

**Problem:** The database file cannot be found.

**Solutions:**
```cpp
// Use absolute path
Thermochimica::setThermoFilename(ctx, "/full/path/to/CO.dat");

// Or ensure working directory is correct
// Run from thermochimica root: ./build/my_app
Thermochimica::setThermoFilename(ctx, "data/CO.dat");
```

### "Temperature out of range" (Error 1)

**Problem:** Temperature is invalid (negative, zero, or NaN).

**Solutions:**
```cpp
// Ensure positive temperature
double T = std::max(1.0, userTemperature);  // At least 1 K
Thermochimica::setTemperature(ctx, T);

// Check for NaN
if (std::isnan(T) || std::isinf(T)) {
    std::cerr << "Invalid temperature\n";
    return;
}
```

### "No species in system" (Error 21)

**Problem:** The specified elements don't match any species in the database.

**Solutions:**
```cpp
// Check what elements the database contains
for (int i = 0; i < Thermochimica::getNumberElementsDatabase(ctx); ++i) {
    std::cout << Thermochimica::getElementAtIndex(ctx, i) << "\n";
}

// Use elements that exist in the database
```

### "GEM solver did not converge" (Error 30)

**Problem:** The solver couldn't find equilibrium within max iterations.

**Solutions:**

1. **Relax tolerances:**
```cpp
ctx.thermo->tolerances[kTolConvergence] = 1e-4;
ctx.thermo->tolerances[kTolDrivingForce] = 1e-3;
```

2. **Try different initial conditions:**
```cpp
// Slight perturbation
Thermochimica::setTemperature(ctx, T + 1.0);
```

3. **Use reinitialization from a nearby solved point:**
```cpp
// Solve at easier conditions first
Thermochimica::setTemperaturePressure(ctx, 1500.0, 1.0);
Thermochimica::thermochimica(ctx);
Thermochimica::saveReinitData(ctx);

// Then move to difficult conditions
Thermochimica::resetThermo(ctx);
Thermochimica::setTemperaturePressure(ctx, difficultT, 1.0);
Thermochimica::setReinitRequested(ctx, true);
Thermochimica::thermochimica(ctx);
```

### "Singular matrix" (Error 32)

**Problem:** Linear algebra solver encountered a singular matrix.

**Solutions:**

1. **Check for zero composition:**
```cpp
// Ensure all specified elements have non-zero mass
if (mass <= 0) {
    std::cerr << "Mass must be positive\n";
}
```

2. **Avoid exact stoichiometric ratios at phase boundaries:**
```cpp
// Instead of exact 1:1
Thermochimica::setElementMass(ctx, 26, 1.0);
Thermochimica::setElementMass(ctx, 8, 1.001);  // Slight offset
```

### "Unsupported phase type" (Error 12)

**Problem:** Database contains a phase model not implemented.

**Solution:** Use a different database or contact maintainers for support.

```cpp
// Check which phase types are in the database
auto& thermo = *ctx.thermo;
for (int i = 0; i < thermo.nSolnPhasesSys; ++i) {
    std::cout << thermo.cSolnPhaseName[i] << ": "
              << thermo.cSolnPhaseType[i] << "\n";
}
```

---

## Debugging Tips

### Enable Verbose Output

```cpp
Thermochimica::setPrintResultsMode(ctx, 2);  // Detailed output
```

### Inspect Internal State

```cpp
if (!ctx.isSuccess()) {
    auto& thermo = *ctx.thermo;
    auto& io = *ctx.io;
    auto& gem = *ctx.gem;

    std::cerr << "Debug info:\n";
    std::cerr << "  Temperature: " << io.dTemperature << " K\n";
    std::cerr << "  Pressure: " << io.dPressure << " atm\n";
    std::cerr << "  nElements: " << thermo.nElements << "\n";
    std::cerr << "  nSpecies: " << thermo.nSpecies << "\n";
    std::cerr << "  Iterations: " << gem.iterGlobal << "\n";
    std::cerr << "  Converged: " << gem.lConverged << "\n";
    std::cerr << "  Function norm: " << gem.dGEMFunctionNorm << "\n";

    std::cerr << "\n  Element masses:\n";
    for (int i = 0; i < 20; ++i) {
        if (io.dElementMass[i] > 0) {
            std::cerr << "    Z=" << i << ": " << io.dElementMass[i] << "\n";
        }
    }
}
```

### Check Solver Progress

```cpp
// Access iteration history
std::cerr << "Solver iterations: " << ctx.gem->iterGlobal << "\n";
std::cerr << "Function norm: " << ctx.gem->dGEMFunctionNorm << "\n";
```

---

## Robust Error Handling Pattern

```cpp
#include <thermochimica/Thermochimica.hpp>
#include <iostream>

bool runCalculation(Thermochimica::ThermoContext& ctx,
                   double T, double P,
                   const std::vector<std::pair<int, double>>& composition) {
    // Reset for new calculation
    Thermochimica::resetThermo(ctx);

    // Validate inputs
    if (T <= 0 || std::isnan(T)) {
        std::cerr << "Invalid temperature: " << T << "\n";
        return false;
    }
    if (P <= 0 || std::isnan(P)) {
        std::cerr << "Invalid pressure: " << P << "\n";
        return false;
    }

    // Set conditions
    Thermochimica::setTemperaturePressure(ctx, T, P);

    for (const auto& [Z, mass] : composition) {
        if (mass > 0) {
            Thermochimica::setElementMass(ctx, Z, mass);
        }
    }

    // Run calculation
    Thermochimica::thermochimica(ctx);

    // Check result
    if (!ctx.isSuccess()) {
        std::cerr << "Calculation failed at T=" << T << " K, P=" << P << " atm\n";
        std::cerr << "Error: " << Thermochimica::getErrorMessage(ctx.infoThermo()) << "\n";
        return false;
    }

    return true;
}
```

---

[← Back to README](../README.md) | [Examples](examples.md) | [Architecture →](architecture.md)
