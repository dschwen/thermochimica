# Architecture

[← Back to README](../README.md) | [Error Handling](error-handling.md) | [API Reference](api-reference.md)

---

Internal design of Thermochimica for contributors and advanced users.

## Design Philosophy

The C++ implementation follows these principles:

1. **Context-based state** - All state is encapsulated in `ThermoContext`, eliminating global variables
2. **Thread safety** - Each context is independent, enabling parallel calculations
3. **RAII memory management** - Smart pointers handle all allocations
4. **Clear separation** - Parser, solver, and output stages are distinct

---

## ThermoContext Structure

The central `ThermoContext` class contains all calculation state:

```cpp
class ThermoContext {
public:
    std::unique_ptr<ThermoState> thermo;   // Core thermodynamic data
    std::unique_ptr<ThermoIO> io;          // Input/output state
    std::unique_ptr<GEMState> gem;         // GEM solver state
    std::unique_ptr<ParserState> parser;   // Parser temporary data
    std::unique_ptr<SubMinState> submin;   // Subminimization state
    std::unique_ptr<CTZState> ctz;         // Common tangent zone
    std::unique_ptr<ReinitState> reinit;   // Reinitialization data
};
```

### ThermoState

Core thermodynamic data (species, phases, equilibrium results):

```cpp
struct ThermoState {
    // Dimensions
    int nElements, nSpecies, nParam;
    int nConPhases, nSolnPhases, nSolnPhasesSys;

    // Species data (Eigen vectors)
    Eigen::VectorXd dStdGibbsEnergy;     // Standard Gibbs energy
    Eigen::VectorXd dMolesSpecies;       // Moles at equilibrium
    Eigen::VectorXd dMolFraction;        // Mole fractions
    Eigen::VectorXd dChemicalPotential;  // Chemical potentials (μ/RT)

    // Phase data
    Eigen::VectorXi iPhase;              // Phase index per species
    Eigen::VectorXi iAssemblage;         // Current phase assemblage
    Eigen::VectorXd dMolesPhase;         // Moles per phase

    // Stoichiometry
    Eigen::MatrixXd dStoichSpecies;      // nSpecies × nElements

    // Names
    std::vector<std::string> cSpeciesName;
    std::vector<std::string> cElementName;
    std::vector<std::string> cSolnPhaseName;

    // Tolerances
    Tolerances tolerances;
};
```

### ThermoIO

Input conditions and output flags:

```cpp
struct ThermoIO {
    // Input
    double dTemperature;                  // Temperature (K, internal)
    double dPressure;                     // Pressure (atm, internal)
    std::array<double, 169> dElementMass; // Element masses by Z
    std::string cThermoFileName;

    // Units
    std::string cInputUnitTemperature;    // "K", "C", "F", "R"
    std::string cInputUnitPressure;       // "atm", "bar", "Pa", etc.
    std::string cInputUnitMass;           // "moles", "grams", etc.

    // Flags
    int iPrintResultsMode;
    bool lWriteJSON;
    bool lHeatCapacityEntropyEnthalpy;

    // Output
    int INFOThermo;                       // Error code
    double dGibbsEnergySys;               // System Gibbs energy
    double dHeatCapacity, dEntropy, dEnthalpy;
};
```

### GEMState

Solver iteration state:

```cpp
struct GEMState {
    int iterGlobal;                // Current iteration
    bool lConverged;               // Convergence flag
    double dGEMFunctionNorm;       // Residual norm
    double dMaxSpeciesChange;      // Max mole fraction change

    Eigen::VectorXd dPartialExcessGibbs;  // Partial molar excess G
    Eigen::VectorXd dUpdateVar;           // Newton update direction
};
```

---

## Calculation Flow

The main `thermochimica()` function executes these stages:

```
┌─────────────────────────────────────────────────────────────┐
│                    thermochimica(ctx)                       │
├─────────────────────────────────────────────────────────────┤
│  1. init(ctx)          - Initialize arrays, reset state    │
│  2. checkSystem(ctx)   - Validate T, P, composition        │
│  3. compThermoData()   - Compute Gibbs energies from T     │
│  4. setup(ctx)         - Leveling solver, initial assembly │
│  5. solve(ctx)         - GEM iteration loop                │
│  6. postProcess(ctx)   - Compute output quantities         │
└─────────────────────────────────────────────────────────────┘
```

### Stage Details

**init()**: Reset solver state, allocate work arrays

**checkSystem()**: Validate inputs
- Temperature > 0
- Pressure > 0
- At least one element with mass > 0
- Unit conversion to internal units (K, atm, moles)

**compThermoData()**: Compute standard Gibbs energies
- Evaluate G(T) polynomials for each species
- Handle temperature range selection

**setup()**: Initial phase assemblage
- Run leveling solver to estimate element potentials
- Determine which phases are likely stable
- Set initial mole fractions

**solve()**: GEM iteration
```
while (not converged and iter < maxIter):
    1. Compute chemical potentials
    2. Check phase stability (driving forces)
    3. Add/remove phases as needed
    4. Compute Newton direction
    5. Line search for step size
    6. Update mole fractions
    7. Check convergence
```

**postProcess()**: Final output
- Compute system Gibbs energy
- Prepare output arrays
- Optionally compute Cp, S, H

---

## GEM Solver Algorithm

The Gibbs Energy Minimization algorithm finds the composition that minimizes:

```
G = Σᵢ nᵢ μᵢ
```

subject to mass balance constraints:

```
Σᵢ aᵢⱼ nᵢ = bⱼ   for each element j
```

### Key Components

**Leveling Solver** (`LevelingSolver.cpp`):
- Estimates initial element potentials
- Uses linear programming approach
- Identifies candidate stable phases

**Newton Direction** (`GEMNewton.cpp`):
- Computes search direction from Hessian
- Handles constrained optimization

**Line Search** (`GEMLineSearch.cpp`):
- Finds optimal step size along Newton direction
- Ensures Gibbs energy decreases

**Phase Assemblage** (`PhaseAssemblage.cpp`):
- Checks driving forces for absent phases
- Adds phases with negative driving force
- Removes phases with vanishing moles

**Subminimization** (`Subminimization.cpp`):
- Optimizes within each solution phase
- Updates site fractions for sublattice models

---

## Phase Models

Solution phase excess Gibbs energy models:

### IDMX (Ideal Mixing)

```
G^ex = 0
G^mix = RT Σᵢ xᵢ ln(xᵢ)
```

No parameters needed. Used for ideal gas mixtures.

### QKTO (Kohler-Toop)

Asymmetric interpolation using chemical groups:

```
G^ex = Σᵢⱼ xᵢ xⱼ Σᵥ Lᵢⱼ^(v) (xᵢ - xⱼ)^v
```

With Kohler interpolation for ternary and higher.

### RKMP (Redlich-Kister-Muggianu)

Polynomial excess model:

```
G^ex = Σᵢⱼ xᵢ xⱼ Σᵥ Lᵢⱼ^(v) (xᵢ - xⱼ)^v
      + Σᵢⱼₖ xᵢ xⱼ xₖ Lᵢⱼₖ + ...
```

### SUBL (Compound Energy Formalism)

Sublattice model for ordered phases:

```
G = Σ_compounds y₁ y₂ ... G°_compound
  + RT Σ_sublattices aₛ Σ_constituents yᵢₛ ln(yᵢₛ)
  + G^ex
```

Site fractions yᵢₛ on each sublattice s.

### SUBG/SUBQ (Modified Quasichemical)

Pair/quadruplet approximation for short-range order:

```
G = Σ_pairs n_ij g°_ij + G^config + G^ex
```

Used for oxide and silicate melts.

---

## File Organization

```
src/
├── api/
│   ├── Thermochimica.cpp      # Public API implementation
│   └── CouplingAPI.cpp        # Coupling interface
├── context/
│   ├── ThermoContext.cpp      # Context management
│   ├── ThermoState.cpp        # State allocation
│   └── ...                    # Other state classes
├── parser/
│   └── ChemSageParser.cpp     # Database parsing
├── setup/
│   ├── InitThermo.cpp         # Initialization
│   └── CheckThermoInput.cpp   # Input validation
├── models/
│   ├── ExcessGibbs.cpp        # Model dispatcher
│   ├── IdealMixing.cpp        # IDMX model
│   ├── QKTO.cpp               # Kohler-Toop
│   ├── RKMP.cpp               # Redlich-Kister
│   ├── SUBL.cpp               # Sublattice CEF
│   └── SUBG.cpp               # Quasichemical
├── solver/
│   ├── GEMSolver.cpp          # Main GEM loop
│   ├── GEMNewton.cpp          # Newton direction
│   ├── GEMLineSearch.cpp      # Line search
│   ├── LevelingSolver.cpp     # Initial estimates
│   ├── PhaseAssemblage.cpp    # Phase add/remove
│   ├── Subminimization.cpp    # Phase optimization
│   └── CheckConvergence.cpp   # Convergence check
└── postprocess/
    ├── PostProcess.cpp        # Output preparation
    ├── HeatCapacity.cpp       # Cp/S/H calculation
    └── PrintResults.cpp       # Text output
```

---

## Adding New Features

### Adding a New Phase Model

1. Create `src/models/NewModel.cpp`
2. Implement `compExcessGibbsEnergyNEW(ctx, phaseIndex)`
3. Add enum value to `Constants::PhaseType`
4. Update `ExcessGibbs.cpp` dispatcher
5. Update `ChemSageParser.cpp` to recognize the type

### Adding a New API Function

1. Declare in `include/thermochimica/Thermochimica.hpp`
2. Implement in `src/api/Thermochimica.cpp`
3. Add tests in `tests/`

### Adding New Tolerances

1. Add enum value to `ToleranceIndex` in `Tolerances.hpp`
2. Set default in `Tolerances::initDefaults()`
3. Update `Constants::kNumTolerances` if needed
4. Use in solver code

---

## Testing

Tests use Google Test framework:

```cpp
TEST(TestGroup, TestName) {
    ThermoContext ctx;
    setStandardUnits(ctx);
    setThermoFilename(ctx, "data/CO.dat");
    parseCSDataFile(ctx);

    setTemperaturePressure(ctx, 1000.0, 1.0);
    setElementMass(ctx, 6, 1.0);

    thermochimica(ctx);

    ASSERT_EQ(ctx.infoThermo(), 0);
    EXPECT_NEAR(getGibbsEnergy(ctx), -629533.0, 1000.0);
}
```

Run tests:
```bash
cd build
ctest --output-on-failure
```

---

[← Back to README](../README.md) | [Error Handling](error-handling.md) | [API Reference](api-reference.md)
