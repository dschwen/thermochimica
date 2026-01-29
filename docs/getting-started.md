# Getting Started

[← Back to README](../README.md) | [API Reference →](api-reference.md)

---

This guide walks you through your first Thermochimica calculation.

## Prerequisites

Ensure you have built the library (see [Building](building.md)):

```bash
cd thermochimica
mkdir build && cd build
cmake ..
make
```

## Your First Calculation

This example calculates the equilibrium composition of C + 2O at 1000 K and 1 atm.

```cpp
#include <thermochimica/Thermochimica.hpp>
#include <iostream>

int main() {
    // 1. Create a context (holds all state)
    Thermochimica::ThermoContext ctx;

    // 2. Set units: Kelvin, atmospheres, moles
    Thermochimica::setStandardUnits(ctx);

    // 3. Load thermodynamic database
    Thermochimica::setThermoFilename(ctx, "CO.dat");
    Thermochimica::parseCSDataFile(ctx);

    if (ctx.infoThermo() != 0) {
        std::cerr << "Failed to load database: "
                  << Thermochimica::getErrorMessage(ctx.infoThermo())
                  << std::endl;
        return 1;
    }

    // 4. Set thermodynamic conditions
    Thermochimica::setTemperaturePressure(ctx, 1000.0, 1.0);

    // 5. Set composition by atomic number
    Thermochimica::setElementMass(ctx, 6, 1.0);   // 1 mol Carbon (Z=6)
    Thermochimica::setElementMass(ctx, 8, 2.0);   // 2 mol Oxygen (Z=8)

    // 6. Run equilibrium calculation
    Thermochimica::thermochimica(ctx);

    // 7. Check results
    if (ctx.isSuccess()) {
        std::cout << "Gibbs Energy: "
                  << Thermochimica::getGibbsEnergy(ctx)
                  << " J" << std::endl;

        Thermochimica::printResults(ctx);
    } else {
        std::cerr << "Calculation failed: "
                  << Thermochimica::getErrorMessage(ctx.infoThermo())
                  << std::endl;
        return 1;
    }

    return 0;
}
```

## Step-by-Step Explanation

### Step 1: Create Context

```cpp
Thermochimica::ThermoContext ctx;
```

The `ThermoContext` holds all calculation state. Each context is independent, enabling thread-safe parallel calculations.

### Step 2: Set Units

```cpp
Thermochimica::setStandardUnits(ctx);  // K, atm, moles
```

Standard units are Kelvin, atmospheres, and moles. For other units:

```cpp
Thermochimica::setUnitTemperature(ctx, "C");    // Celsius
Thermochimica::setUnitPressure(ctx, "bar");     // bar
Thermochimica::setUnitMass(ctx, "grams");       // grams
```

### Step 3: Load Database

```cpp
Thermochimica::setThermoFilename(ctx, "CO.dat");
Thermochimica::parseCSDataFile(ctx);
```

Database files use ChemSage format (`.dat`). See [Databases](databases.md) for available files.

Always check for parse errors:

```cpp
if (ctx.infoThermo() != 0) {
    // Handle error
}
```

### Step 4: Set Conditions

```cpp
Thermochimica::setTemperaturePressure(ctx, 1000.0, 1.0);
```

Temperature and pressure in your chosen units.

### Step 5: Set Composition

```cpp
Thermochimica::setElementMass(ctx, 6, 1.0);   // Carbon by atomic number
Thermochimica::setElementMass(ctx, "O", 2.0); // Oxygen by symbol
```

Elements can be specified by atomic number (1-118) or symbol.

### Step 6: Run Calculation

```cpp
Thermochimica::thermochimica(ctx);
```

This runs the complete GEM solver to find equilibrium.

### Step 7: Get Results

```cpp
if (ctx.isSuccess()) {
    double G = Thermochimica::getGibbsEnergy(ctx);
    // Access more results...
}
```

## Expected Output

For C + 2O at 1000 K, the equilibrium is predominantly CO2 gas:

```
Gibbs Energy: -629533 J

Equilibrium phases:
  gas_ideal: 1.0 mol
    CO2: x = 0.9999
    CO:  x = 0.0001
```

## Next Steps

- [API Reference](api-reference.md) - All available functions
- [Examples](examples.md) - More complex calculations
- [Configuration](configuration.md) - Customize tolerances and options

---

[← Back to README](../README.md) | [API Reference →](api-reference.md)
