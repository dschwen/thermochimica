# API Reference

[← Back to README](../README.md) | [Getting Started](getting-started.md) | [Examples →](examples.md)

---

Complete reference for all public Thermochimica functions.

All functions are in the `Thermochimica` namespace.

## Main Solver Functions

### thermochimica

```cpp
void thermochimica(ThermoContext& ctx);
```

Main entry point. Performs complete equilibrium calculation including initialization, solving, and post-processing.

**Example:**
```cpp
Thermochimica::thermochimica(ctx);
if (ctx.isSuccess()) {
    // Calculation succeeded
}
```

### Individual Solver Steps

For fine-grained control, the calculation can be broken into steps:

```cpp
void init(ThermoContext& ctx);           // Initialize solver
void checkSystem(ThermoContext& ctx);    // Validate inputs
void compThermoData(ThermoContext& ctx); // Compute Gibbs energies
void setup(ThermoContext& ctx);          // Setup phase assemblage
void solve(ThermoContext& ctx);          // Run GEM solver
void postProcess(ThermoContext& ctx);    // Compute output quantities
```

---

## Parser Functions

### parseCSDataFile

```cpp
void parseCSDataFile(ThermoContext& ctx, const std::string& filename);
void parseCSDataFile(ThermoContext& ctx);  // Uses filename from context
```

Parse a ChemSage format database file (`.dat`).

**Example:**
```cpp
Thermochimica::setThermoFilename(ctx, "CO.dat");
Thermochimica::parseCSDataFile(ctx);

// Or directly:
Thermochimica::parseCSDataFile(ctx, "CO.dat");
```

---

## Input Setting Functions

### Temperature and Pressure

```cpp
void setTemperaturePressure(ThermoContext& ctx, double temperature, double pressure);
void setTemperature(ThermoContext& ctx, double temperature);
void setPressure(ThermoContext& ctx, double pressure);
```

Set thermodynamic conditions. Values are interpreted in the current unit system.

**Example:**
```cpp
Thermochimica::setTemperaturePressure(ctx, 1000.0, 1.0);  // 1000 K, 1 atm
```

### Composition

```cpp
void setElementMass(ThermoContext& ctx, int atomicNumber, double mass);
void setElementMass(ThermoContext& ctx, const std::string& elementName, double mass);
```

Set element amounts by atomic number (1-118) or element symbol.

**Example:**
```cpp
Thermochimica::setElementMass(ctx, 6, 1.0);     // 1 mol Carbon
Thermochimica::setElementMass(ctx, "O", 2.0);   // 2 mol Oxygen
Thermochimica::setElementMass(ctx, 26, 0.5);    // 0.5 mol Iron
```

### Units

```cpp
void setUnits(ThermoContext& ctx,
              const std::string& tempUnit,
              const std::string& pressUnit,
              const std::string& massUnit);

void setUnitTemperature(ThermoContext& ctx, const std::string& unit);
void setUnitPressure(ThermoContext& ctx, const std::string& unit);
void setUnitMass(ThermoContext& ctx, const std::string& unit);
void setStandardUnits(ThermoContext& ctx);
```

**Temperature units:** `"K"`, `"C"`, `"F"`, `"R"` (Rankine)

**Pressure units:** `"atm"`, `"psi"`, `"bar"`, `"Pa"`, `"kPa"`

**Mass units:** `"moles"`, `"grams"`, `"kilograms"`, `"kg"`, `"g"`, `"atoms"`, `"atom fraction"`, `"mole fraction"`

**Example:**
```cpp
Thermochimica::setStandardUnits(ctx);  // K, atm, moles

// Or explicitly:
Thermochimica::setUnits(ctx, "C", "bar", "grams");
```

### Filename

```cpp
void setThermoFilename(ThermoContext& ctx, const std::string& filename);
```

Set the database filename before parsing.

### Options

```cpp
void setFuzzyStoich(ThermoContext& ctx, bool enable);
void setFuzzyMagnitude(ThermoContext& ctx, double magnitude);
void setPrintResultsMode(ThermoContext& ctx, int mode);
void setWriteJSON(ThermoContext& ctx, bool enable);
void setHeatCapacityEntropyEnthalpy(ThermoContext& ctx, bool enable);
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

### Setting Constraints

```cpp
void setSolnPhaseConstraint(ThermoContext& ctx,
                            const std::string& phaseName,
                            double targetFraction);

void setCondPhaseConstraint(ThermoContext& ctx,
                            const std::string& speciesName,
                            double targetFraction);

void setPhaseConstraint(ThermoContext& ctx,
                        int phaseIndex,
                        bool isSolutionPhase,
                        double targetFraction);
```

Set a phase fraction constraint. Target fraction must be in range [0, 1].

**Example:**
```cpp
// Constrain FCCN phase to 40% element fraction
Thermochimica::setSolnPhaseConstraint(ctx, "FCCN", 0.4);

// Constrain a pure condensed phase
Thermochimica::setCondPhaseConstraint(ctx, "Graphite", 0.1);

// By index (solution phase)
Thermochimica::setPhaseConstraint(ctx, 1, true, 0.5);
```

### Removing Constraints

```cpp
void removePhaseConstraint(ThermoContext& ctx, const std::string& phaseName);
void clearPhaseConstraints(ThermoContext& ctx);
```

**Example:**
```cpp
Thermochimica::removePhaseConstraint(ctx, "FCCN");  // Remove single constraint
Thermochimica::clearPhaseConstraints(ctx);          // Remove all constraints
```

### Querying Phase Fractions

```cpp
std::pair<double, int> getPhaseElementFraction(const ThermoContext& ctx,
                                               const std::string& phaseName);

bool arePhaseConstraintsSatisfied(const ThermoContext& ctx);
```

**Example:**
```cpp
auto [fraction, info] = Thermochimica::getPhaseElementFraction(ctx, "FCCN");
if (info == 0) {
    std::cout << "FCCN element fraction: " << fraction << "\n";
}

if (Thermochimica::arePhaseConstraintsSatisfied(ctx)) {
    std::cout << "All constraints satisfied\n";
}
```

### Solver Parameters

```cpp
void setConstraintTolerance(ThermoContext& ctx, double tolerance);
void setConstraintPenaltyParameter(ThermoContext& ctx, double rho);
void setConstraintMaxOuterIterations(ThermoContext& ctx, int maxIter);
```

| Function | Default | Description |
|----------|---------|-------------|
| `setConstraintTolerance` | 1e-4 | Tolerance for constraint satisfaction |
| `setConstraintPenaltyParameter` | 1.0 | Initial penalty parameter for augmented Lagrangian |
| `setConstraintMaxOuterIterations` | 20 | Maximum outer iterations |

**Example:**
```cpp
Thermochimica::setConstraintTolerance(ctx, 1e-5);
Thermochimica::setConstraintMaxOuterIterations(ctx, 50);
```

### Complete Example

```cpp
#include <thermochimica/Thermochimica.hpp>

int main() {
    Thermochimica::ThermoContext ctx;

    // Load database
    Thermochimica::setStandardUnits(ctx);
    Thermochimica::setThermoFilename(ctx, "NobleMetals-Kaye.dat");
    Thermochimica::parseCSDataFile(ctx);

    // Set conditions
    Thermochimica::setTemperaturePressure(ctx, 1500.0, 1.0);
    Thermochimica::setElementMass(ctx, 42, 0.5);  // Mo
    Thermochimica::setElementMass(ctx, 44, 0.5);  // Ru

    // Constrain FCCN phase to 40% element fraction
    Thermochimica::setSolnPhaseConstraint(ctx, "FCCN", 0.4);

    // Run constrained equilibrium
    Thermochimica::thermochimica(ctx);

    if (ctx.isSuccess()) {
        auto [frac, info] = Thermochimica::getPhaseElementFraction(ctx, "FCCN");
        std::cout << "FCCN fraction: " << frac << std::endl;
    }

    return 0;
}
```

---

## Output Retrieval Functions

### System Properties

```cpp
double getGibbsEnergy(const ThermoContext& ctx);
double getHeatCapacity(const ThermoContext& ctx);
double getEntropy(const ThermoContext& ctx);
double getEnthalpy(const ThermoContext& ctx);
```

Returns system Gibbs energy (J), heat capacity (J/mol/K), entropy (J/mol/K), and enthalpy (J/mol).

**Note:** Heat capacity, entropy, and enthalpy require `setHeatCapacityEntropyEnthalpy(ctx, true)` before calculation.

### Phase Information

```cpp
std::pair<double, int> getMolesPhase(const ThermoContext& ctx,
                                     const std::string& phaseName);
std::pair<double, int> getSolnPhaseMol(const ThermoContext& ctx,
                                       const std::string& phaseName);
std::pair<double, int> getPureConPhaseMol(const ThermoContext& ctx,
                                          const std::string& phaseName);
int getPhaseIndex(const ThermoContext& ctx, const std::string& phaseName);
bool isPhaseGas(const ThermoContext& ctx, const std::string& phaseName);
bool isPhaseMQM(const ThermoContext& ctx, const std::string& phaseName);
```

**Return format:** `std::pair<value, info>` where `info=0` indicates success.

**Example:**
```cpp
auto [moles, info] = Thermochimica::getMolesPhase(ctx, "gas_ideal");
if (info == 0) {
    std::cout << "Gas phase: " << moles << " mol\n";
}
```

### Species Information

```cpp
std::tuple<double, double, int> getOutputSolnSpecies(
    const ThermoContext& ctx,
    const std::string& phaseName,
    const std::string& speciesName);

std::tuple<double, double, int> getOutputMolSpecies(
    const ThermoContext& ctx,
    const std::string& phaseName,
    const std::string& speciesName);
```

| Function | Returns |
|----------|---------|
| `getOutputSolnSpecies` | `(moleFraction, chemicalPotential, info)` |
| `getOutputMolSpecies` | `(moles, moleFraction, info)` |

**Example:**
```cpp
auto [moleFrac, chemPot, info] = Thermochimica::getOutputSolnSpecies(ctx, "gas_ideal", "CO2");
```

### Chemical Potentials

```cpp
std::pair<double, int> getOutputChemPot(const ThermoContext& ctx,
                                        const std::string& elementName);
```

Get chemical potential of an element.

### Element in Phase

```cpp
std::pair<double, int> getElementMolesInPhase(
    const ThermoContext& ctx,
    const std::string& elementName,
    const std::string& phaseName);
```

### Site Fractions (Sublattice Phases)

```cpp
std::pair<double, int> getOutputSiteFraction(
    const ThermoContext& ctx,
    const std::string& phaseName,
    int sublattice,
    const std::string& constituentName);
```

For SUBL/SUBG/SUBI phase types with sublattice models.

---

## Database Query Functions

```cpp
int getNumberElementsDatabase(const ThermoContext& ctx);
std::string getElementAtIndex(const ThermoContext& ctx, int index);

std::pair<int, int> getNumberPhasesSystem(const ThermoContext& ctx);
std::string getPhaseNameAtIndex(const ThermoContext& ctx, int index);

int getNumberSpeciesSystem(const ThermoContext& ctx);
std::string getSpeciesAtIndex(const ThermoContext& ctx, int index);
```

**Note:** All indices are 0-based.

**Example:**
```cpp
int nElements = Thermochimica::getNumberElementsDatabase(ctx);
for (int i = 0; i < nElements; ++i) {
    std::cout << Thermochimica::getElementAtIndex(ctx, i) << "\n";
}

auto [nSoln, nCond] = Thermochimica::getNumberPhasesSystem(ctx);
std::cout << "Solution phases: " << nSoln << ", Condensed: " << nCond << "\n";
```

---

## Reinitialization Functions

For warm restarts when running sequential calculations:

```cpp
void saveReinitData(ThermoContext& ctx);
void setReinitRequested(ThermoContext& ctx, bool requested);
bool isReinitDataAvailable(const ThermoContext& ctx);

std::tuple<std::vector<int>, std::vector<double>, std::vector<double>,
           std::vector<double>, std::vector<double>>
getReinitData(const ThermoContext& ctx);

void setReinitData(ThermoContext& ctx,
                   const std::vector<int>& assemblage,
                   const std::vector<double>& chemPot,
                   const std::vector<double>& molesPhase,
                   const std::vector<double>& elemPot,
                   const std::vector<double>& molFrac);
```

**Example:**
```cpp
// First calculation
Thermochimica::thermochimica(ctx);
Thermochimica::saveReinitData(ctx);

// Second calculation at nearby conditions
Thermochimica::setTemperaturePressure(ctx, 1005.0, 1.0);
Thermochimica::setReinitRequested(ctx, true);
Thermochimica::thermochimica(ctx);  // Faster convergence
```

---

## Reset Functions

```cpp
void resetThermo(ThermoContext& ctx);     // Reset for new calculation (keeps database)
void resetThermoAll(ThermoContext& ctx);  // Full reset (clears database)
```

**Example:**
```cpp
// Multiple calculations with same database
for (double T = 500; T <= 2000; T += 100) {
    Thermochimica::resetThermo(ctx);  // Clear previous results
    Thermochimica::setTemperaturePressure(ctx, T, 1.0);
    Thermochimica::setElementMass(ctx, 6, 1.0);
    Thermochimica::thermochimica(ctx);
}
```

---

## Output Functions

```cpp
void printResults(ThermoContext& ctx);
void printResultsDetailed(ThermoContext& ctx);
void writeJSON(ThermoContext& ctx, bool append = false);
void computeHeatCapacity(ThermoContext& ctx);
```

---

## Utility Functions

```cpp
int getAtomicNumber(const std::string& symbol);
std::string getElementSymbol(int atomicNumber);
const char* getErrorMessage(int errorCode);
```

**Example:**
```cpp
int Z = Thermochimica::getAtomicNumber("Fe");  // Returns 26
std::string sym = Thermochimica::getElementSymbol(26);  // Returns "Fe"

if (ctx.infoThermo() != 0) {
    std::cerr << Thermochimica::getErrorMessage(ctx.infoThermo());
}
```

---

## ThermoContext Class

```cpp
class ThermoContext {
public:
    std::unique_ptr<ThermoState> thermo;           // Core thermodynamic state
    std::unique_ptr<ThermoIO> io;                  // Input/output state
    std::unique_ptr<GEMState> gem;                 // Solver state
    std::unique_ptr<ParserState> parser;           // Parser state
    std::unique_ptr<SubMinState> submin;           // Subminimization state
    std::unique_ptr<CTZState> ctz;                 // Common tangent zone state
    std::unique_ptr<ReinitState> reinit;           // Reinitialization state
    std::unique_ptr<PhaseConstraints> phaseConstraints;  // Phase fraction constraints

    int infoThermo() const;      // Get error code (0 = success)
    bool isSuccess() const;      // Check if last calculation succeeded
    bool isDatabaseLoaded() const;
};
```

**Direct state access:**
```cpp
// Access temperature
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
