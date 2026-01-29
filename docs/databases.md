# Databases

[← Back to README](../README.md) | [Configuration](configuration.md) | [Examples →](examples.md)

---

Thermochimica reads thermodynamic data from ChemSage format files (`.dat`).

## Available Databases

The `data/` directory contains pre-built databases for various chemical systems:

### Simple Systems

| File | Elements | Description |
|------|----------|-------------|
| `CO.dat` | C, O | Carbon-oxygen gas phase (IDMX) |
| `HO.dat` | H, O | Hydrogen-oxygen system |

### Noble Metals

| File | Elements | Description |
|------|----------|-------------|
| `NobleMetals-Kaye.dat` | Mo, Pd, Ru, Tc | Binary/ternary noble metals |
| `MoPdRuRhTc-Kaye.dat` | Mo, Pd, Ru, Rh, Tc | Extended noble metals |
| `PdRuTcMo.dat` | Pd, Ru, Tc, Mo | Large database with many phases |
| `ternaryMiscibility-Kaye.dat` | Mo, Pd, Ru | Ternary miscibility gap |

### Oxide Systems

| File | Elements | Description |
|------|----------|-------------|
| `FeTiVO.dat` | Fe, Ti, V, O | Iron-titanium-vanadium oxides |
| `WAuArO-1.dat`, `WAuArO-2.dat` | W, Au, Ar, O | Tungsten-gold-argon-oxygen |
| `WAuArNeO-1.dat`, `WAuArNeO-2.dat` | W, Au, Ar, Ne, O | Extended gas systems |

### Metallic Systems

| File | Elements | Description |
|------|----------|-------------|
| `AlMg-Liang.dat` | Al, Mg | Aluminum-magnesium |
| `CuFeC-Kang.dat` | Cu, Fe, C | Copper-iron-carbon |
| `ZrH-Dupin.dat` | Zr, H | Zirconium-hydrogen |

### Ionic/Salt Systems

| File | Elements | Description |
|------|----------|-------------|
| `CsI-Pham.dat` | Cs, I | Cesium iodide (ionic) |
| `CsTe-1.dat`, `CsTe-2.dat` | Cs, Te | Cesium telluride |
| `ClAlNa.dat` | Cl, Al, Na | Chloride-aluminum-sodium |

### Sulfide Systems

| File | Elements | Description |
|------|----------|-------------|
| `CaMnS.dat` | Ca, Mn, S | Calcium-manganese sulfide |
| `FeMnCaS-1.dat`, `FeMnCaS-2.dat` | Fe, Mn, Ca, S | Iron-manganese-calcium sulfide |

### Zirconium Systems

| File | Elements | Description |
|------|----------|-------------|
| `ZIRC_no_liq.dat` | Zr, O, ... | Zirconium (no liquid phase) |
| `ZIRC-noSUBI.dat` | Zr, O, ... | Zirconium (no ionic sublattice) |
| `ZrFeKClNaFOLi.dat` | Zr, Fe, K, Cl, Na, F, O, Li | Complex multi-component |

---

## Phase Model Types

ChemSage databases can contain different solution phase models:

| Type | Code | Description | Example Use |
|------|------|-------------|-------------|
| Ideal Mixing | `IDMX` | Ideal solution (ΔG_mix = RT Σ xᵢ ln xᵢ) | Gas phases |
| Kohler-Toop | `QKTO` | Asymmetric Kohler interpolation | FCC/HCP metals |
| Redlich-Kister | `RKMP` | Polynomial excess Gibbs | Binary alloys |
| RKMP + Magnetic | `RKMPM` | RKMP with magnetic ordering | Fe-Ni alloys |
| Sublattice (CEF) | `SUBL` | Compound Energy Formalism | Intermetallics |
| SUBL + Magnetic | `SUBLM` | CEF with magnetic | Magnetic compounds |
| Quasichemical | `SUBG` | Modified Quasichemical Model | Oxide melts |
| SUBQ variant | `SUBQ` | MQM quadruplet approximation | Silicate melts |
| Ionic Sublattice | `SUBI` | Ionic sublattice model | Molten salts |

### Checking Phase Type

```cpp
bool isGas = Thermochimica::isPhaseGas(ctx, "gas_ideal");
bool isMQM = Thermochimica::isPhaseMQM(ctx, "SLAG");
```

---

## Loading a Database

```cpp
// Set filename
Thermochimica::setThermoFilename(ctx, "data/CO.dat");

// Parse
Thermochimica::parseCSDataFile(ctx);

// Check for errors
if (ctx.infoThermo() != 0) {
    std::cerr << Thermochimica::getErrorMessage(ctx.infoThermo());
    return 1;
}

// Query database contents
int nElem = Thermochimica::getNumberElementsDatabase(ctx);
int nSpec = Thermochimica::getNumberSpeciesSystem(ctx);
auto [nSoln, nCond] = Thermochimica::getNumberPhasesSystem(ctx);

std::cout << "Elements: " << nElem << "\n";
std::cout << "Species: " << nSpec << "\n";
std::cout << "Solution phases: " << nSoln << "\n";
std::cout << "Condensed phases: " << nCond << "\n";
```

---

## ChemSage File Format

ChemSage `.dat` files are structured text files. The general format is:

```
[Header section]
  - Database name and metadata
  - Number of elements, species, phases

[System stoichiometry]
  - Element names and atomic masses
  - Stoichiometry coefficients

[Pure species data]
  - Standard Gibbs energy coefficients
  - Temperature ranges

[Solution phase definitions]
  - Phase names and types
  - Mixing model parameters
  - Excess Gibbs energy coefficients
```

### Gibbs Energy Equations

Species Gibbs energy is typically expressed as:

```
G(T) = a + b*T + c*T*ln(T) + d*T² + e*T³ + f/T + ...
```

Coefficients are stored in the database for each species and temperature range.

---

## Using Custom Databases

To use your own ChemSage database:

1. Ensure it follows ChemSage format
2. Place it in an accessible directory
3. Provide the full or relative path

```cpp
Thermochimica::setThermoFilename(ctx, "/path/to/my/database.dat");
Thermochimica::parseCSDataFile(ctx);
```

### Supported Features

- Multiple temperature ranges per species
- Binary, ternary, and quaternary interaction parameters
- Sublattice phases with multiple sublattices
- Magnetic ordering contributions
- Ionic species and charged constituents

### Limitations

- Some advanced ChemSage features may not be fully supported
- Very large databases may require more memory
- Some specialized phase types may not be implemented

---

## Example: Exploring Database Contents

```cpp
#include <thermochimica/Thermochimica.hpp>
#include <iostream>

int main() {
    Thermochimica::ThermoContext ctx;
    Thermochimica::setStandardUnits(ctx);
    Thermochimica::setThermoFilename(ctx, "data/NobleMetals-Kaye.dat");
    Thermochimica::parseCSDataFile(ctx);

    if (!ctx.isDatabaseLoaded()) {
        std::cerr << "Failed to load database\n";
        return 1;
    }

    // List elements
    std::cout << "Elements:\n";
    for (int i = 0; i < Thermochimica::getNumberElementsDatabase(ctx); ++i) {
        std::cout << "  " << Thermochimica::getElementAtIndex(ctx, i) << "\n";
    }

    // List phases
    auto [nSoln, nCond] = Thermochimica::getNumberPhasesSystem(ctx);
    std::cout << "\nSolution phases (" << nSoln << "):\n";
    for (int i = 0; i < nSoln; ++i) {
        std::cout << "  " << Thermochimica::getPhaseNameAtIndex(ctx, i) << "\n";
    }

    // List first 10 species
    int nSpec = Thermochimica::getNumberSpeciesSystem(ctx);
    std::cout << "\nSpecies (first 10 of " << nSpec << "):\n";
    for (int i = 0; i < std::min(10, nSpec); ++i) {
        std::cout << "  " << Thermochimica::getSpeciesAtIndex(ctx, i) << "\n";
    }

    return 0;
}
```

---

[← Back to README](../README.md) | [Configuration](configuration.md) | [Examples →](examples.md)
