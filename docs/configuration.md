# Configuration

[← Back to README](../README.md) | [Building](building.md) | [Databases →](databases.md)

---

Customize tolerances, units, and solver options.

## Tolerances

Thermochimica uses 15 solver tolerances that can be adjusted for specific calculations.

### Accessing Tolerances

```cpp
#include <thermochimica/util/Tolerances.hpp>

// Set a tolerance
ctx.thermo->tolerances[kTolMassBalance] = 1e-6;

// Get a tolerance
double tol = ctx.thermo->tolerances[kTolConvergence];
```

### Tolerance Reference

| Index | Name | Default | Description |
|-------|------|---------|-------------|
| 0 | `kTolMassBalance` | 1e-5 | Mass balance residual tolerance |
| 1 | `kTolChemPotential` | 1e-6 | Chemical potential residual |
| 2 | `kTolGibbsEnergy` | 1e-8 | Gibbs energy tolerance |
| 3 | `kTolMoleFraction` | 1e-12 | Minimum mole fraction |
| 4 | `kTolDrivingForce` | 1e-4 | Driving force for phase stability |
| 5 | `kTolElementPotential` | 1e-6 | Element potential tolerance |
| 6 | `kTolConvergence` | 1e-6 | Overall convergence criterion |
| 7 | `kTolPhaseChange` | 1e-3 | Phase addition/removal threshold |
| 8 | `kTolSiteFraction` | 1e-6 | Site fraction tolerance (sublattice) |
| 9 | `kTolMinMoles` | 1e-20 | Minimum moles threshold |
| 10 | `kTolPhaseMoles` | 1e-9 | Phase moles tolerance |
| 11 | `kTolFunctionNorm` | 1e-6 | Function norm tolerance |
| 12 | `kTolRelative` | 1e-8 | Relative tolerance |
| 13 | `kTolAbsolute` | 1e-10 | Absolute tolerance |
| 14 | `kTolEuclidean` | 1e-2 | Euclidean norm tolerance |

### Example: Tightening Convergence

```cpp
// Tighter convergence for high-precision calculations
ctx.thermo->tolerances[kTolConvergence] = 1e-8;
ctx.thermo->tolerances[kTolMassBalance] = 1e-7;
ctx.thermo->tolerances[kTolChemPotential] = 1e-8;
```

### Example: Relaxing for Difficult Systems

```cpp
// More relaxed tolerances for difficult convergence
ctx.thermo->tolerances[kTolConvergence] = 1e-4;
ctx.thermo->tolerances[kTolDrivingForce] = 1e-3;
```

### Reset to Defaults

```cpp
ctx.thermo->tolerances.initDefaults();
```

---

## Units

### Temperature Units

| Unit | String | Description |
|------|--------|-------------|
| Kelvin | `"K"` | SI unit (default) |
| Celsius | `"C"` | T(K) = T(C) + 273.15 |
| Fahrenheit | `"F"` | T(K) = (T(F) + 459.67) × 5/9 |
| Rankine | `"R"` | T(K) = T(R) × 5/9 |

```cpp
Thermochimica::setUnitTemperature(ctx, "C");
Thermochimica::setTemperature(ctx, 726.85);  // = 1000 K
```

### Pressure Units

| Unit | String | Description |
|------|--------|-------------|
| Atmosphere | `"atm"` | Default |
| Bar | `"bar"` | 1 bar = 0.986923 atm |
| Pascal | `"Pa"` | 1 Pa = 9.86923e-6 atm |
| Kilopascal | `"kPa"` | 1 kPa = 0.00986923 atm |
| PSI | `"psi"` | 1 psi = 0.068046 atm |

```cpp
Thermochimica::setUnitPressure(ctx, "bar");
Thermochimica::setPressure(ctx, 1.01325);  // = 1 atm
```

### Mass Units

| Unit | Strings | Description |
|------|---------|-------------|
| Moles | `"moles"`, `"mol"` | Default |
| Grams | `"grams"`, `"g"` | Converted using atomic mass |
| Kilograms | `"kilograms"`, `"kg"` | 1000 grams |
| Atoms | `"atoms"` | Converted via Avogadro |
| Atom fraction | `"atom fraction"` | Normalized to 1 |
| Mole fraction | `"mole fraction"` | Normalized to 1 |

```cpp
Thermochimica::setUnitMass(ctx, "grams");
Thermochimica::setElementMass(ctx, 26, 55.845);  // 55.845 g Fe = 1 mol
```

### Standard Units

```cpp
Thermochimica::setStandardUnits(ctx);  // K, atm, moles
```

### Combined Unit Setting

```cpp
Thermochimica::setUnits(ctx, "C", "bar", "grams");
```

---

## Solver Options

### Print Results Mode

Control output verbosity:

```cpp
Thermochimica::setPrintResultsMode(ctx, 0);  // No output
Thermochimica::setPrintResultsMode(ctx, 1);  // Basic output
Thermochimica::setPrintResultsMode(ctx, 2);  // Detailed output
```

### Heat Capacity, Entropy, Enthalpy

Enable calculation of thermodynamic properties:

```cpp
Thermochimica::setHeatCapacityEntropyEnthalpy(ctx, true);
Thermochimica::thermochimica(ctx);

double Cp = Thermochimica::getHeatCapacity(ctx);
double S = Thermochimica::getEntropy(ctx);
double H = Thermochimica::getEnthalpy(ctx);
```

### JSON Output

Enable JSON file output:

```cpp
Thermochimica::setWriteJSON(ctx, true);
Thermochimica::thermochimica(ctx);
// Creates thermoout.json
```

### Fuzzy Stoichiometry

For approximate stoichiometry matching:

```cpp
Thermochimica::setFuzzyStoich(ctx, true);
Thermochimica::setFuzzyMagnitude(ctx, 0.01);  // 1% tolerance
```

---

## Direct State Access

For advanced usage, you can access internal state directly.

### ThermoIO (Input/Output State)

```cpp
// Input conditions
ctx.io->dTemperature;           // Temperature (K)
ctx.io->dPressure;              // Pressure (atm)
ctx.io->dElementMass[Z];        // Element mass by atomic number

// Input settings
ctx.io->cThermoFileName;        // Database filename
ctx.io->iPrintResultsMode;      // Print mode (0, 1, 2)
ctx.io->lWriteJSON;             // JSON output flag
ctx.io->lHeatCapacityEntropyEnthalpy;  // Enable Cp/S/H

// Output
ctx.io->INFOThermo;             // Error code
ctx.io->dGibbsEnergySys;        // System Gibbs energy
ctx.io->dHeatCapacity;          // Heat capacity
ctx.io->dEntropy;               // Entropy
ctx.io->dEnthalpy;              // Enthalpy
```

### ThermoState (Core State)

```cpp
// System dimensions
ctx.thermo->nElements;          // Number of elements
ctx.thermo->nSpecies;           // Number of species
ctx.thermo->nSolnPhasesSys;     // Number of solution phases
ctx.thermo->nConPhasesSys;      // Number of condensed phases

// Arrays (Eigen vectors/matrices)
ctx.thermo->dStdGibbsEnergy(i);    // Gibbs energy of species i
ctx.thermo->dMolFraction(i);       // Mole fraction of species i
ctx.thermo->dMolesSpecies(i);      // Moles of species i
ctx.thermo->dChemicalPotential(i); // Chemical potential
ctx.thermo->dStoichSpecies(i, j);  // Stoichiometry matrix

// Names
ctx.thermo->cSpeciesName[i];       // Species name
ctx.thermo->cElementName[j];       // Element name
ctx.thermo->cSolnPhaseName[k];     // Phase name
```

### GEMState (Solver State)

```cpp
ctx.gem->iterGlobal;            // Current iteration count
ctx.gem->lConverged;            // Convergence flag
ctx.gem->dGEMFunctionNorm;      // Function norm
```

---

## Example: Custom Configuration

```cpp
#include <thermochimica/Thermochimica.hpp>

int main() {
    Thermochimica::ThermoContext ctx;

    // Custom units
    Thermochimica::setUnits(ctx, "C", "bar", "grams");

    // Tighter tolerances
    ctx.thermo->tolerances[kTolConvergence] = 1e-8;
    ctx.thermo->tolerances[kTolMassBalance] = 1e-7;

    // Enable all output
    Thermochimica::setPrintResultsMode(ctx, 2);
    Thermochimica::setHeatCapacityEntropyEnthalpy(ctx, true);
    Thermochimica::setWriteJSON(ctx, true);

    // Load and run
    Thermochimica::setThermoFilename(ctx, "CO.dat");
    Thermochimica::parseCSDataFile(ctx);

    Thermochimica::setTemperaturePressure(ctx, 726.85, 1.01325);  // 1000 K, 1 atm
    Thermochimica::setElementMass(ctx, 6, 12.011);   // 12.011 g C
    Thermochimica::setElementMass(ctx, 8, 32.0);     // 32.0 g O

    Thermochimica::thermochimica(ctx);

    return 0;
}
```

---

[← Back to README](../README.md) | [Building](building.md) | [Databases →](databases.md)
