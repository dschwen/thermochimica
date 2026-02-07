# Examples

[← Back to README](../README.md) | [Databases](databases.md) | [Error Handling →](error-handling.md)

---

Code examples for common Thermochimica use cases using the modern `ThermoClass` API.

## Basic Equilibrium Calculation

Calculate equilibrium for C + 2O at 1000 K:

```cpp
#include <thermochimica/ThermoClass.hpp>
#include <iostream>

int main() {
    using namespace Thermochimica;

    ThermoClass thermo;

    // Setup
    thermo.setStandardUnits();
    thermo.loadDatabase("data/CO.dat");

    // Conditions
    thermo.setTemperaturePressure(1000.0, 1.0);
    thermo.setElementMass(6, 1.0);   // C
    thermo.setElementMass(8, 2.0);   // O

    // Calculate
    thermo.calculate();

    // Results
    if (thermo.isSuccess()) {
        std::cout << "Gibbs Energy: " << thermo.getGibbsEnergy() << " J\n";
        thermo.printResults();
    }

    return 0;
}
```

---

## Multi-Point Calculations (Temperature Sweep)

Calculate equilibrium at multiple temperatures:

```cpp
#include <thermochimica/ThermoClass.hpp>
#include <iostream>
#include <iomanip>

int main() {
    using namespace Thermochimica;

    ThermoClass thermo;

    thermo.setStandardUnits();
    thermo.loadDatabase("data/CO.dat");

    std::cout << std::setw(10) << "T (K)" << std::setw(20) << "G (J)\n";
    std::cout << std::string(30, '-') << "\n";

    for (double T = 500; T <= 2000; T += 100) {
        // Reset for new calculation (keeps database loaded)
        thermo.reset();

        thermo.setTemperaturePressure(T, 1.0);
        thermo.setElementMass(6, 1.0);
        thermo.setElementMass(8, 2.0);

        thermo.calculate();

        if (thermo.isSuccess()) {
            std::cout << std::setw(10) << T
                      << std::setw(20) << thermo.getGibbsEnergy()
                      << "\n";
        }
    }

    return 0;
}
```

---

## Accessing Internal State

Direct access to species mole fractions and chemical potentials:

```cpp
#include <thermochimica/ThermoClass.hpp>
#include <iostream>
#include <iomanip>

int main() {
    using namespace Thermochimica;

    ThermoClass thermo;

    thermo.setStandardUnits();
    thermo.loadDatabase("data/CO.dat");

    thermo.setTemperaturePressure(1000.0, 1.0);
    thermo.setElementMass(6, 1.0);
    thermo.setElementMass(8, 2.0);

    thermo.calculate();

    if (thermo.isSuccess()) {
        // Access internal state via context
        ThermoContext& ctx = thermo.getContext();
        auto& state = *ctx.thermo;

        std::cout << "Species results:\n";
        std::cout << std::setw(15) << "Species"
                  << std::setw(15) << "Mole Frac"
                  << std::setw(15) << "Moles"
                  << std::setw(15) << "mu/RT\n";
        std::cout << std::string(60, '-') << "\n";

        for (int i = 0; i < state.nSpecies; ++i) {
            if (state.dMolFraction(i) > 1e-10) {
                std::cout << std::setw(15) << state.cSpeciesName[i]
                          << std::setw(15) << state.dMolFraction(i)
                          << std::setw(15) << state.dMolesSpecies(i)
                          << std::setw(15) << state.dChemicalPotential(i)
                          << "\n";
            }
        }

        std::cout << "\nElement potentials:\n";
        for (int j = 0; j < state.nElements; ++j) {
            std::cout << "  " << state.cElementName[j] << ": "
                      << state.dElementPotential(j) << "\n";
        }
    }

    return 0;
}
```

---

## Warm Restart (Reinitialization)

Use previous results as initial guess for faster convergence:

```cpp
#include <thermochimica/ThermoClass.hpp>
#include <iostream>

int main() {
    using namespace Thermochimica;

    ThermoClass thermo;

    thermo.setStandardUnits();
    thermo.loadDatabase("data/NobleMetals-Kaye.dat");

    // First calculation
    thermo.setTemperaturePressure(1500.0, 1.0);
    thermo.setElementMass(42, 0.5);  // Mo
    thermo.setElementMass(44, 0.5);  // Ru

    thermo.calculate();

    if (thermo.isSuccess()) {
        std::cout << "T=1500 K: G = " << thermo.getGibbsEnergy() << " J\n";

        // Save state for reinitialization
        thermo.saveReinitData();
    }

    // Second calculation at nearby temperature (faster)
    thermo.reset();
    thermo.setTemperaturePressure(1510.0, 1.0);
    thermo.setElementMass(42, 0.5);
    thermo.setElementMass(44, 0.5);

    // Request reinitialization from saved data
    thermo.setReinitRequested(true);
    thermo.calculate();

    if (thermo.isSuccess()) {
        std::cout << "T=1510 K: G = " << thermo.getGibbsEnergy() << " J\n";
    }

    return 0;
}
```

---

## Thread-Safe Parallel Calculations

Each ThermoClass instance is independent, enabling safe parallel execution:

```cpp
#include <thermochimica/ThermoClass.hpp>
#include <iostream>
#include <thread>
#include <vector>
#include <mutex>

std::mutex outputMutex;

void calculateAtTemperature(double T, const std::string& dataFile) {
    using namespace Thermochimica;

    // Each thread has its own instance
    ThermoClass thermo;

    thermo.setStandardUnits();
    thermo.loadDatabase(dataFile);

    thermo.setTemperaturePressure(T, 1.0);
    thermo.setElementMass(6, 1.0);
    thermo.setElementMass(8, 2.0);

    thermo.calculate();

    if (thermo.isSuccess()) {
        std::lock_guard<std::mutex> lock(outputMutex);
        std::cout << "T=" << T << " K: G = "
                  << thermo.getGibbsEnergy() << " J\n";
    }
}

int main() {
    std::vector<std::thread> threads;

    // Launch parallel calculations
    for (double T = 500; T <= 2000; T += 250) {
        threads.emplace_back(calculateAtTemperature, T, "data/CO.dat");
    }

    // Wait for all threads
    for (auto& t : threads) {
        t.join();
    }

    return 0;
}
```

---

## Heat Capacity, Entropy, and Enthalpy

Calculate thermodynamic properties:

```cpp
#include <thermochimica/ThermoClass.hpp>
#include <iostream>
#include <iomanip>

int main() {
    using namespace Thermochimica;

    ThermoClass thermo;

    thermo.setStandardUnits();
    thermo.loadDatabase("data/CO.dat");

    // Enable property calculation
    thermo.setHeatCapacityEntropyEnthalpy(true);

    thermo.setTemperaturePressure(1000.0, 1.0);
    thermo.setElementMass(6, 1.0);
    thermo.setElementMass(8, 2.0);

    thermo.calculate();

    if (thermo.isSuccess()) {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "Gibbs Energy:  " << thermo.getGibbsEnergy() << " J\n";
        std::cout << "Heat Capacity: " << thermo.getHeatCapacity() << " J/(mol K)\n";
        std::cout << "Entropy:       " << thermo.getEntropy() << " J/(mol K)\n";
        std::cout << "Enthalpy:      " << thermo.getEnthalpy() << " J/mol\n";
    }

    return 0;
}
```

---

## Custom Tolerances

Adjust solver tolerances for specific needs:

```cpp
#include <thermochimica/ThermoClass.hpp>
#include <thermochimica/util/Tolerances.hpp>
#include <iostream>

int main() {
    using namespace Thermochimica;

    ThermoClass thermo;

    thermo.setStandardUnits();
    thermo.loadDatabase("data/CO.dat");

    // Access context for tolerance adjustment
    ThermoContext& ctx = thermo.getContext();

    // Tighten convergence for high precision
    ctx.thermo->tolerances[kTolConvergence] = 1e-10;
    ctx.thermo->tolerances[kTolMassBalance] = 1e-8;
    ctx.thermo->tolerances[kTolChemPotential] = 1e-10;

    thermo.setTemperaturePressure(1000.0, 1.0);
    thermo.setElementMass(6, 1.0);
    thermo.setElementMass(8, 2.0);

    thermo.calculate();

    if (thermo.isSuccess()) {
        std::cout << "Converged with tight tolerances\n";
        std::cout << "Gibbs Energy: " << std::setprecision(12)
                  << thermo.getGibbsEnergy() << " J\n";
    }

    return 0;
}
```

---

## Using Different Units

Work with Celsius and bar instead of Kelvin and atm:

```cpp
#include <thermochimica/ThermoClass.hpp>
#include <iostream>

int main() {
    using namespace Thermochimica;

    ThermoClass thermo;

    // Set custom units
    thermo.setUnitTemperature("C");
    thermo.setUnitPressure("bar");
    thermo.setUnitMass("grams");

    thermo.loadDatabase("data/CO.dat");

    // Now use Celsius and bar
    thermo.setTemperaturePressure(726.85, 1.01325);  // ~1000 K, ~1 atm

    // Mass in grams
    thermo.setElementMass(6, 12.011);   // 12.011 g C = 1 mol
    thermo.setElementMass(8, 32.0);     // 32.0 g O = 2 mol

    thermo.calculate();

    if (thermo.isSuccess()) {
        std::cout << "Gibbs Energy: " << thermo.getGibbsEnergy() << " J\n";
    }

    return 0;
}
```

---

## Phase and Species Queries

Get detailed information about equilibrium phases:

```cpp
#include <thermochimica/ThermoClass.hpp>
#include <iostream>
#include <iomanip>

int main() {
    using namespace Thermochimica;

    ThermoClass thermo;

    thermo.setStandardUnits();
    thermo.loadDatabase("data/NobleMetals-Kaye.dat");

    thermo.setTemperaturePressure(1800.0, 1.0);
    thermo.setElementMass(42, 0.33);  // Mo
    thermo.setElementMass(44, 0.33);  // Ru
    thermo.setElementMass(46, 0.34);  // Pd

    thermo.calculate();

    if (thermo.isSuccess()) {
        auto [nSoln, nCond] = thermo.getNumberPhasesSystem();

        std::cout << "Equilibrium phases:\n\n";

        // Check each phase
        for (int i = 0; i < nSoln; ++i) {
            std::string phaseName = thermo.getPhaseNameAtIndex(i);
            auto [moles, info] = thermo.getMolesPhase(phaseName);

            if (info == 0 && moles > 1e-10) {
                std::cout << phaseName << ": " << moles << " mol\n";

                // Get species in this phase
                for (int j = 0; j < thermo.getNumberSpeciesSystem(); ++j) {
                    std::string specName = thermo.getSpeciesAtIndex(j);
                    auto [molFrac, chemPot, sinfo] =
                        thermo.getOutputSolnSpecies(phaseName, specName);

                    if (sinfo == 0 && molFrac > 1e-6) {
                        std::cout << "  " << specName << ": x = " << molFrac << "\n";
                    }
                }
                std::cout << "\n";
            }
        }
    }

    return 0;
}
```

---

## Phase Fraction Constraints

Constrain phase fractions for phase field modeling applications.

### Fixed Assemblage Mode

When constraints are active, Thermochimica uses a **fixed assemblage mode**:

- **Constrained phases are forced into the assemblage** - regardless of thermodynamic stability
- **Unconstrained phases are excluded** - not included in the calculation
- **Constraints must sum to ~1.0** - all element mass must be accounted for

This is designed for phase field coupling where the mesoscale simulation imposes phase fractions, and Thermochimica computes the resulting chemical potentials and compositions.

### Two-Phase Constraint Example

```cpp
#include <thermochimica/ThermoClass.hpp>
#include <iostream>

int main() {
    using namespace Thermochimica;

    ThermoClass thermo;

    thermo.setStandardUnits();
    thermo.loadDatabase("data/NobleMetals-Kaye.dat");

    // Set conditions
    thermo.setTemperaturePressure(1500.0, 1.0);
    thermo.setElementMass(42, 0.5);  // Mo
    thermo.setElementMass(44, 0.5);  // Ru

    // Constrain TWO phases with fractions summing to 1.0
    // This is required: all element mass must go into constrained phases
    thermo.setSolnPhaseConstraint("FCCN", 0.4);
    thermo.setSolnPhaseConstraint("BCCN", 0.6);

    // Optionally adjust solver parameters
    thermo.setConstraintTolerance(1e-4);
    thermo.setConstraintMaxOuterIterations(30);

    thermo.calculate();

    if (thermo.isSuccess()) {
        auto [fcc_frac, fcc_info] = thermo.getPhaseElementFraction("FCCN");
        auto [bcc_frac, bcc_info] = thermo.getPhaseElementFraction("BCCN");
        std::cout << "FCCN fraction: " << fcc_frac << "\n";  // ~0.4
        std::cout << "BCCN fraction: " << bcc_frac << "\n";  // ~0.6

        // Check constraint satisfaction
        if (thermo.arePhaseConstraintsSatisfied()) {
            std::cout << "All phase constraints satisfied.\n";
        }
    }

    // Clean up constraints for next calculation
    thermo.clearPhaseConstraints();

    return 0;
}
```

### Important Notes

**Constraints must account for all mass:**
```cpp
// CORRECT: fractions sum to 1.0
thermo.setSolnPhaseConstraint("FCCN", 0.4);
thermo.setSolnPhaseConstraint("BCCN", 0.6);

// INCORRECT: single constraint without accounting for remaining mass
// The phase will end up at 1.0 (100%) since it's the only phase present
thermo.setSolnPhaseConstraint("FCCN", 0.4);  // Will actually be 1.0!
```

**Managing constraints:**
```cpp
// Remove a single constraint
thermo.removePhaseConstraint("BCCN");

// Clear all constraints (returns to unconstrained equilibrium mode)
thermo.clearPhaseConstraints();
```

**Phase fraction definition:** Phase fraction is defined at the element level as `(sum of element moles in phase) / (total element moles in system)`. This differs from mole fraction, making it suitable for phase field applications where you track the spatial distribution of phases.

---

## Custom Solvers (Advanced)

Use the strategy pattern to inject custom solver implementations:

```cpp
#include <thermochimica/ThermoClass.hpp>
#include <thermochimica/interfaces/ISolver.hpp>
#include <iostream>

// Custom solver implementation
class MyCustomSolver : public ISolver {
public:
    int solve(ThermoState& state, ThermoIO& io, GEMState& gemState,
             PhaseAssemblageManager& phaseManager, INewtonSolver& newton,
             ILineSearch& lineSearch,
             const std::vector<IThermodynamicModel*>& models) override {
        // Custom solver logic here...
        std::cout << "Using custom solver!\n";
        // For this example, just return success
        return 0;
    }

    void initialize(ThermoState& state, GEMState& gemState) override {
        std::cout << "Custom solver initialization\n";
    }

    bool isConverged(const ThermoState& state,
                    const GEMState& gemState) const override {
        return gemState.lConverged;
    }

    const char* getSolverName() const override {
        return "MyCustomSolver";
    }
};

int main() {
    using namespace Thermochimica;

    ThermoClass thermo;

    thermo.loadDatabase("data/CO.dat");
    thermo.setStandardUnits();

    // Inject custom solver
    thermo.setSolver(std::make_unique<MyCustomSolver>());

    thermo.setTemperaturePressure(1000.0, 1.0);
    thermo.setElementMass(6, 1.0);
    thermo.setElementMass(8, 2.0);

    thermo.calculate();  // Uses custom solver

    return 0;
}
```

---

## Move Semantics (Modern C++)

Efficient resource management with move semantics:

```cpp
#include <thermochimica/ThermoClass.hpp>
#include <iostream>
#include <vector>

using namespace Thermochimica;

// Factory function returning ThermoClass by value (efficient move)
ThermoClass createConfiguredThermo(const std::string& database) {
    ThermoClass thermo;
    thermo.loadDatabase(database);
    thermo.setStandardUnits();
    return thermo;  // Efficient move, not copy
}

int main() {
    // Move construction
    ThermoClass thermo = createConfiguredThermo("data/CO.dat");

    // Use configured instance
    thermo.setTemperaturePressure(1000.0, 1.0);
    thermo.setElementMass(6, 1.0);
    thermo.setElementMass(8, 2.0);
    thermo.calculate();

    // Store in containers (also uses move)
    std::vector<ThermoClass> thermos;
    thermos.push_back(std::move(thermo));  // Move into vector

    return 0;
}
```

---

[← Back to README](../README.md) | [Databases](databases.md) | [Error Handling →](error-handling.md)
