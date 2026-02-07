# API Reference

[← Back to README](../README.md) | [Getting Started](getting-started.md) | [Examples →](examples.md)

---

Complete reference for the modern `ThermoClass` API.

All classes and functions are in the `Thermochimica` namespace.

**Note:** For migrating from the old free-function API, see [Migration Guide](MIGRATION_FROM_FREE_FUNCTIONS.md).

## ThermoClass Overview

The `ThermoClass` is the modern, canonical API for Thermochimica. All calculations are performed through member functions of this class.

```cpp
#include <thermochimica/ThermoClass.hpp>

Thermochimica::ThermoClass thermo;
```

## Main Solver Functions

### calculate

```cpp
int calculate();
```

Main entry point. Performs complete equilibrium calculation including initialization, solving, and post-processing.

**Returns:** Error code (0 = success)

**Example:**
```cpp
thermo.calculate();
if (thermo.isSuccess()) {
    // Calculation succeeded
}
```

### Individual Solver Steps

For fine-grained control, the calculation can be broken into steps:

```cpp
void initialize();       // Initialize solver
void checkSystem();      // Validate inputs
void computeThermoData();  // Compute Gibbs energies
void setup();            // Setup phase assemblage
int solve();             // Run GEM solver
```

---

## Database Loading

### loadDatabase

```cpp
int loadDatabase(const std::string& filename);
```

Load and parse a ChemSage format database file (`.dat`). Combines `setThermoFilename()` and `parseCSDataFile()`.

**Returns:** Error code (0 = success)

**Example:**
```cpp
int result = thermo.loadDatabase("CO.dat");
if (result != 0) {
    std::cerr << "Failed: " << thermo.getErrorMessage() << "\n";
}
```

### Individual Steps

```cpp
void setThermoFilename(const std::string& filename);
int parseCSDataFile();
```

For cases where you need separate control over filename setting and parsing.

---

## Input Configuration

### Temperature and Pressure

```cpp
void setTemperaturePressure(double temperature, double pressure);
void setTemperature(double temperature);
void setPressure(double pressure);
```

Set thermodynamic conditions. Values are interpreted in the current unit system.

**Example:**
```cpp
thermo.setTemperaturePressure(1000.0, 1.0);  // 1000 K, 1 atm
```

### Composition

```cpp
void setElementMass(int atomicNumber, double mass);
void setElementMass(const std::string& elementName, double mass);
```

Set element amounts by atomic number (1-118) or element symbol.

**Example:**
```cpp
thermo.setElementMass(6, 1.0);     // 1 mol Carbon
thermo.setElementMass("O", 2.0);   // 2 mol Oxygen
thermo.setElementMass(26, 0.5);    // 0.5 mol Iron
```

### Units

```cpp
void setStandardUnits();
void setUnitsSI();
void setUnitTemperature(const std::string& unit);
void setUnitPressure(const std::string& unit);
void setUnitMass(const std::string& unit);
```

**Temperature units:** `"K"`, `"C"`, `"F"`, `"R"` (Rankine)

**Pressure units:** `"atm"`, `"psi"`, `"bar"`, `"Pa"`, `"kPa"`

**Mass units:** `"moles"`, `"grams"`, `"kilograms"`, `"kg"`, `"g"`, `"atoms"`, `"atom fraction"`, `"mole fraction"`

**Example:**
```cpp
thermo.setStandardUnits();  // K, atm, moles

// Or explicitly:
thermo.setUnitTemperature("C");
thermo.setUnitPressure("bar");
thermo.setUnitMass("grams");
```

### Options

```cpp
void setFuzzyStoich(bool enable);
void setFuzzyMagnitude(double magnitude);
void setPrintResultsMode(int mode);
void setWriteJSON(bool enable);
void setHeatCapacityEntropyEnthalpy(bool enable);
```

| Function | Description |
|----------|-------------|
| `setFuzzyStoich` | Enable fuzzy stoichiometry matching |
| `setFuzzyMagnitude` | Set fuzzy matching tolerance |
| `setPrintResultsMode` | 0=none, 1=basic, 2=detailed |
| `setWriteJSON` | Enable JSON output |
| `setHeatCapacityEntropyEnthalpy` | Enable Cp, S, H calculation |

---

## Phase Constraint Functions

For phase field modeling applications, you can constrain phase fractions at the element level.
The phase fraction is defined as: `f_p = (sum of element moles in phase p) / (total element moles in system)`.

### Fixed Assemblage Mode

When constraints are active, Thermochimica uses a **fixed assemblage mode** designed for phase field coupling:

- **Constrained phases are forced into the assemblage** - The solver includes exactly the phases you constrain, regardless of whether they would be thermodynamically stable
- **Unconstrained phases are excluded** - Phases without constraints are not included in the assemblage
- **Constraints should sum to ~1.0** - Since all element mass must go somewhere, constrain all phases that should be present with fractions that sum to 1.0

This approach differs from the classical unconstrained equilibrium calculation where the solver determines which phases are stable. In phase field applications, the phase fractions are imposed by the mesoscale simulation, and Thermochimica computes the corresponding chemical potentials and compositions.

**Important**: Constraint enforcement only works for phases in the assemblage. If you constrain a single phase to 0.4 without constraining other phases, that phase will end up at 1.0 (100%) because it's the only phase present to hold the element mass.

### Setting Constraints

```cpp
void setSolnPhaseConstraint(const std::string& phaseName, double targetFraction);
void setCondPhaseConstraint(const std::string& speciesName, double targetFraction);
```

Set a phase fraction constraint. Target fraction must be in range [0, 1].

**Example:**
```cpp
// Constrain FCCN phase to 40% element fraction
thermo.setSolnPhaseConstraint("FCCN", 0.4);

// Constrain a pure condensed phase
thermo.setCondPhaseConstraint("Graphite", 0.1);
```

### Removing Constraints

```cpp
void removePhaseConstraint(const std::string& phaseName);
void clearPhaseConstraints();
```

**Example:**
```cpp
thermo.removePhaseConstraint("FCCN");  // Remove single constraint
thermo.clearPhaseConstraints();         // Remove all constraints
```

### Querying Phase Fractions

```cpp
std::pair<double, int> getPhaseElementFraction(const std::string& phaseName) const;
bool arePhaseConstraintsSatisfied() const;
```

**Example:**
```cpp
auto [fraction, info] = thermo.getPhaseElementFraction("FCCN");
if (info == 0) {
    std::cout << "FCCN element fraction: " << fraction << "\n";
}

if (thermo.arePhaseConstraintsSatisfied()) {
    std::cout << "All constraints satisfied\n";
}
```

### Solver Parameters

```cpp
void setConstraintTolerance(double tolerance);
void setConstraintPenaltyParameter(double rho);
void setConstraintMaxOuterIterations(int maxIter);
```

| Function | Default | Description |
|----------|---------|-------------|
| `setConstraintTolerance` | 1e-4 | Tolerance for constraint satisfaction |
| `setConstraintPenaltyParameter` | 1.0 | Initial penalty parameter for augmented Lagrangian |
| `setConstraintMaxOuterIterations` | 20 | Maximum outer iterations |

**Example:**
```cpp
thermo.setConstraintTolerance(1e-5);
thermo.setConstraintMaxOuterIterations(50);
```

### Complete Example

```cpp
#include <thermochimica/ThermoClass.hpp>

int main() {
    using namespace Thermochimica;

    ThermoClass thermo;

    // Load database
    thermo.setStandardUnits();
    thermo.loadDatabase("NobleMetals-Kaye.dat");

    // Set conditions
    thermo.setTemperaturePressure(1500.0, 1.0);
    thermo.setElementMass(42, 0.5);  // Mo
    thermo.setElementMass(44, 0.5);  // Ru

    // Constrain two phases with fractions summing to 1.0
    // (Required: all element mass must be accounted for)
    thermo.setSolnPhaseConstraint("FCCN", 0.4);
    thermo.setSolnPhaseConstraint("BCCN", 0.6);

    // Run constrained equilibrium
    thermo.calculate();

    if (thermo.isSuccess()) {
        auto [fcc_frac, fcc_info] = thermo.getPhaseElementFraction("FCCN");
        auto [bcc_frac, bcc_info] = thermo.getPhaseElementFraction("BCCN");
        std::cout << "FCCN fraction: " << fcc_frac << std::endl;  // ~0.4
        std::cout << "BCCN fraction: " << bcc_frac << std::endl;  // ~0.6
    }

    return 0;
}
```

---

## Output Retrieval Functions

### System Properties

```cpp
double getGibbsEnergy() const;
double getHeatCapacity() const;
double getEntropy() const;
double getEnthalpy() const;
```

Returns system Gibbs energy (J), heat capacity (J/mol/K), entropy (J/mol/K), and enthalpy (J/mol).

**Note:** Heat capacity, entropy, and enthalpy require `setHeatCapacityEntropyEnthalpy(true)` before calculation.

### Phase Information

```cpp
std::pair<double, int> getMolesPhase(const std::string& phaseName) const;
std::pair<double, int> getSolnPhaseMol(const std::string& phaseName) const;
std::pair<double, int> getPureConPhaseMol(const std::string& phaseName) const;
int getPhaseIndex(const std::string& phaseName) const;
bool isPhaseGas(const std::string& phaseName) const;
bool isPhaseMQM(const std::string& phaseName) const;
```

**Return format:** `std::pair<value, info>` where `info=0` indicates success.

**Example:**
```cpp
auto [moles, info] = thermo.getMolesPhase("gas_ideal");
if (info == 0) {
    std::cout << "Gas phase: " << moles << " mol\n";
}
```

### Species Information

```cpp
std::tuple<double, double, int> getOutputSolnSpecies(
    const std::string& phaseName,
    const std::string& speciesName) const;

std::tuple<double, double, int> getOutputMolSpecies(
    const std::string& phaseName,
    const std::string& speciesName) const;
```

| Function | Returns |
|----------|---------|
| `getOutputSolnSpecies` | `(moleFraction, chemicalPotential, info)` |
| `getOutputMolSpecies` | `(moles, moleFraction, info)` |

**Example:**
```cpp
auto [moleFrac, chemPot, info] = thermo.getOutputSolnSpecies("gas_ideal", "CO2");
```

### Chemical Potentials

```cpp
std::pair<double, int> getOutputChemPot(const std::string& elementName) const;
```

Get chemical potential of an element.

### Element in Phase

```cpp
std::pair<double, int> getElementMolesInPhase(
    const std::string& elementName,
    const std::string& phaseName) const;
```

### Site Fractions (Sublattice Phases)

```cpp
std::pair<double, int> getOutputSiteFraction(
    const std::string& phaseName,
    int sublattice,
    const std::string& constituentName) const;
```

For SUBL/SUBG/SUBI phase types with sublattice models.

---

## Database Query Functions

```cpp
int getNumberElementsDatabase() const;
std::string getElementAtIndex(int index) const;

std::pair<int, int> getNumberPhasesSystem() const;
std::string getPhaseNameAtIndex(int index) const;

int getNumberSpeciesSystem() const;
std::string getSpeciesAtIndex(int index) const;
```

**Note:** All indices are 0-based.

**Example:**
```cpp
int nElements = thermo.getNumberElementsDatabase();
for (int i = 0; i < nElements; ++i) {
    std::cout << thermo.getElementAtIndex(i) << "\n";
}

auto [nSoln, nCond] = thermo.getNumberPhasesSystem();
std::cout << "Solution phases: " << nSoln << ", Condensed: " << nCond << "\n";
```

---

## Reinitialization Functions

For warm restarts when running sequential calculations:

```cpp
void saveReinitData();
void setReinitRequested(bool requested);
bool isReinitDataAvailable() const;
```

**Example:**
```cpp
// First calculation
thermo.calculate();
thermo.saveReinitData();

// Second calculation at nearby conditions
thermo.reset();
thermo.setTemperaturePressure(1005.0, 1.0);
thermo.setElementMass(6, 1.0);
thermo.setReinitRequested(true);
thermo.calculate();  // Faster convergence
```

---

## Reset Functions

```cpp
void reset();      // Reset for new calculation (keeps database)
void resetAll();   // Full reset (clears database)
```

**Example:**
```cpp
// Multiple calculations with same database
for (double T = 500; T <= 2000; T += 100) {
    thermo.reset();  // Clear previous results
    thermo.setTemperaturePressure(T, 1.0);
    thermo.setElementMass(6, 1.0);
    thermo.calculate();
}
```

---

## Output Functions

```cpp
void printResults() const;
void printResultsDetailed() const;
void writeJSON(bool append = false) const;
void computeHeatCapacity();
```

---

## Status and Error Handling

```cpp
int getInfoCode() const;       // Get error/info code (0 = success)
bool isSuccess() const;        // Check if calculation succeeded
std::string getErrorMessage() const;  // Get error description
```

**Example:**
```cpp
thermo.calculate();

if (!thermo.isSuccess()) {
    std::cerr << "Error code: " << thermo.getInfoCode() << "\n";
    std::cerr << "Error: " << thermo.getErrorMessage() << "\n";
}
```

## Utility Functions (Static)

```cpp
static int getAtomicNumber(const std::string& symbol);
static std::string getElementSymbol(int atomicNumber);
```

**Example:**
```cpp
int Z = ThermoClass::getAtomicNumber("Fe");  // Returns 26
std::string sym = ThermoClass::getElementSymbol(26);  // Returns "Fe"
```

---

## Advanced: Custom Solvers

ThermoClass supports dependency injection for custom implementations:

```cpp
void setSolver(std::unique_ptr<ISolver> solver);
void setNewtonSolver(std::unique_ptr<INewtonSolver> newton);
void setLineSearch(std::unique_ptr<ILineSearch> lineSearch);
```

See [class_based_api.md](class_based_api.md) for examples of custom solver implementations.

## Advanced: Direct State Access

For advanced use cases requiring direct state manipulation:

```cpp
ThermoContext& getContext();
const ThermoContext& getContext() const;
```

**Example:**
```cpp
// Access underlying context
ThermoContext& ctx = thermo.getContext();

// Direct state access
double T = ctx.io->dTemperature;

// Access species data
for (int i = 0; i < ctx.thermo->nSpecies; ++i) {
    std::string name = ctx.thermo->cSpeciesName[i];
    double G = ctx.thermo->dStdGibbsEnergy(i);
    double x = ctx.thermo->dMolFraction(i);
}
```

---

[← Back to README](../README.md) | [Getting Started](getting-started.md) | [Examples →](examples.md)
