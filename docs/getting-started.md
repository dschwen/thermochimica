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
#include <thermochimica/ThermoClass.hpp>
#include <iostream>

int main() {
    using namespace Thermochimica;

    // 1. Create Thermochimica instance
    ThermoClass thermo;

    // 2. Set units: Kelvin, atmospheres, moles
    thermo.setStandardUnits();

    // 3. Load thermodynamic database
    int result = thermo.loadDatabase("CO.dat");

    if (result != 0) {
        std::cerr << "Failed to load database: "
                  << thermo.getErrorMessage()
                  << std::endl;
        return 1;
    }

    // 4. Set thermodynamic conditions
    thermo.setTemperaturePressure(1000.0, 1.0);

    // 5. Set composition by atomic number
    thermo.setElementMass(6, 1.0);   // 1 mol Carbon (Z=6)
    thermo.setElementMass(8, 2.0);   // 2 mol Oxygen (Z=8)

    // 6. Run equilibrium calculation
    thermo.calculate();

    // 7. Check results
    if (thermo.isSuccess()) {
        std::cout << "Gibbs Energy: "
                  << thermo.getGibbsEnergy()
                  << " J" << std::endl;

        thermo.printResults();
    } else {
        std::cerr << "Calculation failed: "
                  << thermo.getErrorMessage()
                  << std::endl;
        return 1;
    }

    return 0;
}
```

## Step-by-Step Explanation

### Step 1: Create ThermoClass Instance

```cpp
Thermochimica::ThermoClass thermo;
```

The `ThermoClass` encapsulates all calculation state. Each instance is independent, enabling thread-safe parallel calculations.

### Step 2: Set Units

```cpp
thermo.setStandardUnits();  // K, atm, moles
```

Standard units are Kelvin, atmospheres, and moles. For other units:

```cpp
thermo.setUnitTemperature("C");    // Celsius
thermo.setUnitPressure("bar");     // bar
thermo.setUnitMass("grams");       // grams
```

### Step 3: Load Database

```cpp
int result = thermo.loadDatabase("CO.dat");
```

Database files use ChemSage format (`.dat`). See [Databases](databases.md) for available files.

Always check for parse errors:

```cpp
if (result != 0) {
    // Handle error
}
```

### Step 4: Set Conditions

```cpp
thermo.setTemperaturePressure(1000.0, 1.0);
```

Temperature and pressure in your chosen units.

### Step 5: Set Composition

```cpp
thermo.setElementMass(6, 1.0);   // Carbon by atomic number
thermo.setElementMass("O", 2.0); // Oxygen by symbol
```

Elements can be specified by atomic number (1-118) or symbol.

### Step 6: Run Calculation

```cpp
thermo.calculate();
```

This runs the complete GEM solver to find equilibrium.

### Step 7: Get Results

```cpp
if (thermo.isSuccess()) {
    double G = thermo.getGibbsEnergy();
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
