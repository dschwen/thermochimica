# ThermoClass API Documentation

## Overview

Thermochimica provides a modern, object-oriented C++ API through the `ThermoClass` class. This is the canonical and recommended interface for all applications.

## Key Benefits

- **Object-Oriented Design**: Natural C++ class-based interface
- **Encapsulation**: All state contained within `ThermoClass` instance
- **Thread-Safe**: Multiple independent instances can run in parallel
- **Extensibility**: Strategy pattern allows custom solvers, models, and line search algorithms
- **Type Safety**: Member functions with compile-time checking
- **RAII Semantics**: Automatic resource management
- **Move Semantics**: Efficient resource transfer

## Basic Usage

### Simple Calculation

```cpp
#include <thermochimica/ThermoClass.hpp>
#include <iostream>

int main() {
    using namespace Thermochimica;

    // Create Thermochimica instance
    ThermoClass thermo;

    // Load thermodynamic database
    thermo.loadDatabase("CO.dat");

    // Set units and conditions
    thermo.setStandardUnits();  // K, atm, moles
    thermo.setTemperaturePressure(1000.0, 1.0);

    // Set composition (atomic number, mass)
    thermo.setElementMass(6, 1.0);   // Carbon: 1.0 moles
    thermo.setElementMass(8, 2.0);   // Oxygen: 2.0 moles

    // Run calculation
    int result = thermo.calculate();

    if (result == 0) {
        // Get results
        double gibbs = thermo.getGibbsEnergy();
        auto [moles, info] = thermo.getMolesPhase("gas_ideal");
        auto [chemPot, info2] = thermo.getOutputChemPot("C");

        std::cout << "Calculation succeeded!" << std::endl;
        std::cout << "Gibbs energy: " << gibbs << " J" << std::endl;
        std::cout << "Gas phase moles: " << moles << std::endl;
        thermo.printResults();
    } else {
        std::cerr << "Error: " << thermo.getErrorMessage() << std::endl;
    }

    return 0;
}
```

### Using Element Names

```cpp
ThermoClass thermo;
thermo.loadDatabase("CO.dat");
thermo.setStandardUnits();
thermo.setTemperaturePressure(1500.0, 1.0);

// Set composition by element name
thermo.setElementMass("C", 1.0);
thermo.setElementMass("O", 1.0);

thermo.calculate();
```

### SI Units

```cpp
ThermoClass thermo;
thermo.loadDatabase("CO.dat");

// Use SI units: K, Pa, moles
thermo.setUnitsSI();

thermo.setTemperature(1000.0);  // K
thermo.setPressure(101325.0);    // Pa (1 atm)
thermo.setElementMass("C", 1.0);
thermo.setElementMass("O", 2.0);

thermo.calculate();
```

## Phase Constraints

Thermochimica supports constraining phase fractions for phase-field modeling applications.

### Solution Phase Constraints

```cpp
ThermoClass thermo;
thermo.loadDatabase("NobleMetals-Kaye.dat");
thermo.setStandardUnits();
thermo.setTemperaturePressure(2000.0, 1.0);
thermo.setElementMass("Mo", 0.5);
thermo.setElementMass("Ru", 0.5);

// Constrain HCP_A3 phase to 30% of elements
thermo.setSolnPhaseConstraint("HCP_A3", 0.3);

thermo.calculate();

// Remove constraint for next calculation
thermo.removePhaseConstraint("HCP_A3");

// Or clear all constraints
thermo.clearPhaseConstraints();
```

### Condensed Phase Constraints

```cpp
ThermoClass thermo;
thermo.loadDatabase("CO.dat");
thermo.setStandardUnits();
thermo.setTemperaturePressure(1000.0, 1.0);
thermo.setElementMass("C", 1.0);
thermo.setElementMass("O", 2.0);

// Constrain graphite to 50% of carbon
thermo.setCondPhaseConstraint("C(gr)", 0.5);

thermo.calculate();
```

## Multiple Independent Calculations

```cpp
// Thread-safe: each instance is independent
ThermoClass thermo1, thermo2;

// Configure first instance
thermo1.loadDatabase("CO.dat");
thermo1.setStandardUnits();
thermo1.setTemperaturePressure(1000.0, 1.0);
thermo1.setElementMass("C", 1.0);
thermo1.setElementMass("O", 2.0);

// Configure second instance with different conditions
thermo2.loadDatabase("CO.dat");
thermo2.setStandardUnits();
thermo2.setTemperaturePressure(1500.0, 2.0);
thermo2.setElementMass("C", 2.0);
thermo2.setElementMass("O", 1.0);

// Run both calculations
int result1 = thermo1.calculate();
int result2 = thermo2.calculate();

// Results are independent
double gibbs1 = thermo1.getGibbsEnergy();
double gibbs2 = thermo2.getGibbsEnergy();
```

## Reusing Instances

```cpp
ThermoClass thermo;
thermo.loadDatabase("CO.dat");
thermo.setStandardUnits();

// First calculation
thermo.setTemperaturePressure(1000.0, 1.0);
thermo.setElementMass("C", 1.0);
thermo.setElementMass("O", 2.0);
thermo.calculate();
double gibbs1 = thermo.getGibbsEnergy();

// Reset for new calculation (keeps database)
thermo.reset();

// Second calculation with different conditions
thermo.setTemperaturePressure(1500.0, 1.0);
thermo.setElementMass("C", 2.0);
thermo.setElementMass("O", 1.0);
thermo.calculate();
double gibbs2 = thermo.getGibbsEnergy();

// Full reset (clears database too)
thermo.resetAll();
```

## Error Handling

```cpp
ThermoClass thermo;

// Check status after each operation
int result = thermo.loadDatabase("nonexistent.dat");
if (result != 0) {
    std::cerr << "Load failed: " << thermo.getErrorMessage() << std::endl;
    return 1;
}

// Alternative: check isSuccess()
thermo.calculate();
if (!thermo.isSuccess()) {
    std::cerr << "Calculation failed: " << thermo.getErrorMessage() << std::endl;
    std::cerr << "Error code: " << thermo.getInfoCode() << std::endl;
    return 1;
}
```

## Advanced: Granular Control

For advanced users who need fine-grained control over the calculation process:

```cpp
ThermoClass thermo;
thermo.loadDatabase("CO.dat");
thermo.setStandardUnits();
thermo.setTemperaturePressure(1000.0, 1.0);
thermo.setElementMass("C", 1.0);
thermo.setElementMass("O", 2.0);

// Manual step-by-step execution
thermo.initialize();       // Initialize solver state
thermo.checkSystem();      // Validate input
thermo.computeThermoData();  // Compute Gibbs energies
thermo.setup();            // Prepare for iteration
int result = thermo.solve();  // Run GEM solver

// Or use the all-in-one calculate() method
// thermo.calculate();
```

## Advanced: Custom Solvers (Extensibility)

The architecture supports custom solver strategies through dependency injection:

```cpp
#include <thermochimica/interfaces/ISolver.hpp>
#include <thermochimica/solver/StandardGEMSolver.hpp>

// Custom solver implementation
class MyCustomSolver : public ISolver {
public:
    int solve(ThermoState& state, ThermoIO& io, GEMState& gemState,
             PhaseAssemblageManager& phaseManager, INewtonSolver& newton,
             ILineSearch& lineSearch,
             const std::vector<IThermodynamicModel*>& models) override {
        // Custom solver logic...
        return 0;
    }

    void initialize(ThermoState& state, GEMState& gemState) override {
        // Custom initialization...
    }

    bool isConverged(const ThermoState& state,
                    const GEMState& gemState) const override {
        // Custom convergence check...
        return gemState.lConverged;
    }

    const char* getSolverName() const override {
        return "MyCustomSolver";
    }
};

// Use custom solver
ThermoClass thermo;
thermo.loadDatabase("CO.dat");

// Inject custom solver
thermo.setSolver(std::make_unique<MyCustomSolver>());

// Or inject custom line search
thermo.setLineSearch(std::make_unique<MyCustomLineSearch>());

thermo.calculate();  // Uses custom solver
```

## Advanced: Accessing Internal State

For advanced use cases requiring direct state access, use `getContext()`:

```cpp
ThermoClass thermo;
thermo.loadDatabase("CO.dat");
thermo.setStandardUnits();
thermo.setTemperaturePressure(1000.0, 1.0);

// Access underlying context for direct state manipulation
ThermoContext& ctx = thermo.getContext();

// Direct access to internal state
double T = ctx.io->dTemperature;
auto& state = *ctx.thermo;
for (int i = 0; i < state.nSpecies; ++i) {
    double x = state.dMolFraction(i);
}
```

## API Reference

### Database Loading
- `void setThermoFilename(const std::string& filename)`
- `int parseCSDataFile()`
- `int loadDatabase(const std::string& filename)` - Combines setThermoFilename + parseCSDataFile

### Input Configuration
- `void setTemperature(double temperature)`
- `void setPressure(double pressure)`
- `void setTemperaturePressure(double temperature, double pressure)`
- `void setElementMass(int atomicNumber, double mass)`
- `void setElementMass(const std::string& elementName, double mass)`
- `void setStandardUnits()` - K, atm, moles
- `void setUnitsSI()` - K, Pa, moles

### Phase Constraints
- `void setSolnPhaseConstraint(const std::string& phaseName, double targetFraction)`
- `void setCondPhaseConstraint(const std::string& speciesName, double targetFraction)`
- `void removePhaseConstraint(const std::string& phaseName)`
- `void clearPhaseConstraints()`

### Calculation
- `int calculate()` - All-in-one calculation
- `void initialize()` - Initialize solver (granular control)
- `void checkSystem()` - Validate input (granular control)
- `void computeThermoData()` - Compute Gibbs energies (granular control)
- `void setup()` - Prepare for iteration (granular control)
- `int solve()` - Run GEM solver (granular control)

### Output Retrieval
- `std::pair<double, int> getOutputChemPot(const std::string& elementName) const`
- `std::pair<double, int> getMolesPhase(const std::string& phaseName) const`
- `double getGibbsEnergy() const`
- `void printResults() const`

### Status
- `int getInfoCode() const` - Get error/info code (0 = success)
- `bool isSuccess() const` - Check if calculation succeeded
- `std::string getErrorMessage() const` - Get error description

### Reset
- `void reset()` - Reset for new calculation (keeps database)
- `void resetAll()` - Reset everything including database

### Advanced Configuration
- `void setSolver(std::unique_ptr<ISolver> solver)` - Inject custom solver
- `void setNewtonSolver(std::unique_ptr<INewtonSolver> newton)` - Inject custom Newton solver
- `void setLineSearch(std::unique_ptr<ILineSearch> lineSearch)` - Inject custom line search

### Legacy Compatibility
- `ThermoContext& getContext()` - Get underlying context for legacy code
- `const ThermoContext& getContext() const` - Get const context

## Migration from Old Code

If migrating from older code, see the comprehensive [Migration Guide](MIGRATION_FROM_FREE_FUNCTIONS.md) for:
- Complete API mapping tables
- Common migration patterns
- Before/after examples
- Troubleshooting tips

## Performance

The class-based API has negligible overhead compared to the free-function API:
- Virtual function call overhead: < 1% (one virtual call per strategy component)
- Memory overhead: Minimal (stores strategy pointers)
- **All 90 tests produce identical results to free-function API**

## Thread Safety

Each `ThermoClass` instance is independent and thread-safe when used from separate threads. However, a single instance should not be accessed from multiple threads simultaneously without external synchronization.

```cpp
// Thread-safe usage
std::vector<std::thread> threads;
for (int i = 0; i < 4; ++i) {
    threads.emplace_back([i]() {
        ThermoClass thermo;  // Each thread gets its own instance
        thermo.loadDatabase("CO.dat");
        thermo.setStandardUnits();
        thermo.setTemperaturePressure(1000.0 + i * 100, 1.0);
        thermo.setElementMass("C", 1.0);
        thermo.setElementMass("O", 2.0);
        thermo.calculate();
    });
}

for (auto& t : threads) {
    t.join();
}
```

## See Also

- [Architecture Documentation](oo_architecture.md) - Design patterns and architecture
- [Examples](../examples/) - Complete example programs
- [API Reference](api_reference.md) - Detailed API documentation
