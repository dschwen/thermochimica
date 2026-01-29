# Examples

[← Back to README](../README.md) | [Databases](databases.md) | [Error Handling →](error-handling.md)

---

Code examples for common Thermochimica use cases.

## Basic Equilibrium Calculation

Calculate equilibrium for C + 2O at 1000 K:

```cpp
#include <thermochimica/Thermochimica.hpp>
#include <iostream>

int main() {
    Thermochimica::ThermoContext ctx;

    // Setup
    Thermochimica::setStandardUnits(ctx);
    Thermochimica::setThermoFilename(ctx, "data/CO.dat");
    Thermochimica::parseCSDataFile(ctx);

    // Conditions
    Thermochimica::setTemperaturePressure(ctx, 1000.0, 1.0);
    Thermochimica::setElementMass(ctx, 6, 1.0);   // C
    Thermochimica::setElementMass(ctx, 8, 2.0);   // O

    // Calculate
    Thermochimica::thermochimica(ctx);

    // Results
    if (ctx.isSuccess()) {
        std::cout << "Gibbs Energy: " << Thermochimica::getGibbsEnergy(ctx) << " J\n";
        Thermochimica::printResults(ctx);
    }

    return 0;
}
```

---

## Multi-Point Calculations (Temperature Sweep)

Calculate equilibrium at multiple temperatures:

```cpp
#include <thermochimica/Thermochimica.hpp>
#include <iostream>
#include <iomanip>

int main() {
    Thermochimica::ThermoContext ctx;

    Thermochimica::setStandardUnits(ctx);
    Thermochimica::setThermoFilename(ctx, "data/CO.dat");
    Thermochimica::parseCSDataFile(ctx);

    std::cout << std::setw(10) << "T (K)" << std::setw(20) << "G (J)\n";
    std::cout << std::string(30, '-') << "\n";

    for (double T = 500; T <= 2000; T += 100) {
        // Reset for new calculation (keeps database loaded)
        Thermochimica::resetThermo(ctx);

        Thermochimica::setTemperaturePressure(ctx, T, 1.0);
        Thermochimica::setElementMass(ctx, 6, 1.0);
        Thermochimica::setElementMass(ctx, 8, 2.0);

        Thermochimica::thermochimica(ctx);

        if (ctx.isSuccess()) {
            std::cout << std::setw(10) << T
                      << std::setw(20) << Thermochimica::getGibbsEnergy(ctx)
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
#include <thermochimica/Thermochimica.hpp>
#include <iostream>
#include <iomanip>

int main() {
    Thermochimica::ThermoContext ctx;

    Thermochimica::setStandardUnits(ctx);
    Thermochimica::setThermoFilename(ctx, "data/CO.dat");
    Thermochimica::parseCSDataFile(ctx);

    Thermochimica::setTemperaturePressure(ctx, 1000.0, 1.0);
    Thermochimica::setElementMass(ctx, 6, 1.0);
    Thermochimica::setElementMass(ctx, 8, 2.0);

    Thermochimica::thermochimica(ctx);

    if (ctx.isSuccess()) {
        auto& thermo = *ctx.thermo;

        std::cout << "Species results:\n";
        std::cout << std::setw(15) << "Species"
                  << std::setw(15) << "Mole Frac"
                  << std::setw(15) << "Moles"
                  << std::setw(15) << "mu/RT\n";
        std::cout << std::string(60, '-') << "\n";

        for (int i = 0; i < thermo.nSpecies; ++i) {
            if (thermo.dMolFraction(i) > 1e-10) {
                std::cout << std::setw(15) << thermo.cSpeciesName[i]
                          << std::setw(15) << thermo.dMolFraction(i)
                          << std::setw(15) << thermo.dMolesSpecies(i)
                          << std::setw(15) << thermo.dChemicalPotential(i)
                          << "\n";
            }
        }

        std::cout << "\nElement potentials:\n";
        for (int j = 0; j < thermo.nElements; ++j) {
            std::cout << "  " << thermo.cElementName[j] << ": "
                      << thermo.dElementPotential(j) << "\n";
        }
    }

    return 0;
}
```

---

## Warm Restart (Reinitialization)

Use previous results as initial guess for faster convergence:

```cpp
#include <thermochimica/Thermochimica.hpp>
#include <iostream>

int main() {
    Thermochimica::ThermoContext ctx;

    Thermochimica::setStandardUnits(ctx);
    Thermochimica::setThermoFilename(ctx, "data/NobleMetals-Kaye.dat");
    Thermochimica::parseCSDataFile(ctx);

    // First calculation
    Thermochimica::setTemperaturePressure(ctx, 1500.0, 1.0);
    Thermochimica::setElementMass(ctx, 42, 0.5);  // Mo
    Thermochimica::setElementMass(ctx, 44, 0.5);  // Ru

    Thermochimica::thermochimica(ctx);

    if (ctx.isSuccess()) {
        std::cout << "T=1500 K: G = " << Thermochimica::getGibbsEnergy(ctx) << " J\n";

        // Save state for reinitialization
        Thermochimica::saveReinitData(ctx);
    }

    // Second calculation at nearby temperature (faster)
    Thermochimica::resetThermo(ctx);
    Thermochimica::setTemperaturePressure(ctx, 1510.0, 1.0);
    Thermochimica::setElementMass(ctx, 42, 0.5);
    Thermochimica::setElementMass(ctx, 44, 0.5);

    // Request reinitialization from saved data
    Thermochimica::setReinitRequested(ctx, true);
    Thermochimica::thermochimica(ctx);

    if (ctx.isSuccess()) {
        std::cout << "T=1510 K: G = " << Thermochimica::getGibbsEnergy(ctx) << " J\n";
    }

    return 0;
}
```

---

## Thread-Safe Parallel Calculations

Each ThermoContext is independent, enabling safe parallel execution:

```cpp
#include <thermochimica/Thermochimica.hpp>
#include <iostream>
#include <thread>
#include <vector>
#include <mutex>

std::mutex outputMutex;

void calculateAtTemperature(double T, const std::string& dataFile) {
    // Each thread has its own context
    Thermochimica::ThermoContext ctx;

    Thermochimica::setStandardUnits(ctx);
    Thermochimica::setThermoFilename(ctx, dataFile);
    Thermochimica::parseCSDataFile(ctx);

    Thermochimica::setTemperaturePressure(ctx, T, 1.0);
    Thermochimica::setElementMass(ctx, 6, 1.0);
    Thermochimica::setElementMass(ctx, 8, 2.0);

    Thermochimica::thermochimica(ctx);

    if (ctx.isSuccess()) {
        std::lock_guard<std::mutex> lock(outputMutex);
        std::cout << "T=" << T << " K: G = "
                  << Thermochimica::getGibbsEnergy(ctx) << " J\n";
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
#include <thermochimica/Thermochimica.hpp>
#include <iostream>
#include <iomanip>

int main() {
    Thermochimica::ThermoContext ctx;

    Thermochimica::setStandardUnits(ctx);
    Thermochimica::setThermoFilename(ctx, "data/CO.dat");
    Thermochimica::parseCSDataFile(ctx);

    // Enable property calculation
    Thermochimica::setHeatCapacityEntropyEnthalpy(ctx, true);

    Thermochimica::setTemperaturePressure(ctx, 1000.0, 1.0);
    Thermochimica::setElementMass(ctx, 6, 1.0);
    Thermochimica::setElementMass(ctx, 8, 2.0);

    Thermochimica::thermochimica(ctx);

    if (ctx.isSuccess()) {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "Gibbs Energy:  " << Thermochimica::getGibbsEnergy(ctx) << " J\n";
        std::cout << "Heat Capacity: " << Thermochimica::getHeatCapacity(ctx) << " J/(mol K)\n";
        std::cout << "Entropy:       " << Thermochimica::getEntropy(ctx) << " J/(mol K)\n";
        std::cout << "Enthalpy:      " << Thermochimica::getEnthalpy(ctx) << " J/mol\n";
    }

    return 0;
}
```

---

## Custom Tolerances

Adjust solver tolerances for specific needs:

```cpp
#include <thermochimica/Thermochimica.hpp>
#include <thermochimica/util/Tolerances.hpp>
#include <iostream>

int main() {
    Thermochimica::ThermoContext ctx;

    Thermochimica::setStandardUnits(ctx);
    Thermochimica::setThermoFilename(ctx, "data/CO.dat");
    Thermochimica::parseCSDataFile(ctx);

    // Tighten convergence for high precision
    ctx.thermo->tolerances[Thermochimica::kTolConvergence] = 1e-10;
    ctx.thermo->tolerances[Thermochimica::kTolMassBalance] = 1e-8;
    ctx.thermo->tolerances[Thermochimica::kTolChemPotential] = 1e-10;

    Thermochimica::setTemperaturePressure(ctx, 1000.0, 1.0);
    Thermochimica::setElementMass(ctx, 6, 1.0);
    Thermochimica::setElementMass(ctx, 8, 2.0);

    Thermochimica::thermochimica(ctx);

    if (ctx.isSuccess()) {
        std::cout << "Converged with tight tolerances\n";
        std::cout << "Gibbs Energy: " << std::setprecision(12)
                  << Thermochimica::getGibbsEnergy(ctx) << " J\n";
    }

    return 0;
}
```

---

## Using Different Units

Work with Celsius and bar instead of Kelvin and atm:

```cpp
#include <thermochimica/Thermochimica.hpp>
#include <iostream>

int main() {
    Thermochimica::ThermoContext ctx;

    // Set custom units
    Thermochimica::setUnitTemperature(ctx, "C");
    Thermochimica::setUnitPressure(ctx, "bar");
    Thermochimica::setUnitMass(ctx, "grams");

    Thermochimica::setThermoFilename(ctx, "data/CO.dat");
    Thermochimica::parseCSDataFile(ctx);

    // Now use Celsius and bar
    Thermochimica::setTemperaturePressure(ctx, 726.85, 1.01325);  // ~1000 K, ~1 atm

    // Mass in grams
    Thermochimica::setElementMass(ctx, 6, 12.011);   // 12.011 g C = 1 mol
    Thermochimica::setElementMass(ctx, 8, 32.0);     // 32.0 g O = 2 mol

    Thermochimica::thermochimica(ctx);

    if (ctx.isSuccess()) {
        std::cout << "Gibbs Energy: " << Thermochimica::getGibbsEnergy(ctx) << " J\n";
    }

    return 0;
}
```

---

## Phase and Species Queries

Get detailed information about equilibrium phases:

```cpp
#include <thermochimica/Thermochimica.hpp>
#include <iostream>
#include <iomanip>

int main() {
    Thermochimica::ThermoContext ctx;

    Thermochimica::setStandardUnits(ctx);
    Thermochimica::setThermoFilename(ctx, "data/NobleMetals-Kaye.dat");
    Thermochimica::parseCSDataFile(ctx);

    Thermochimica::setTemperaturePressure(ctx, 1800.0, 1.0);
    Thermochimica::setElementMass(ctx, 42, 0.33);  // Mo
    Thermochimica::setElementMass(ctx, 44, 0.33);  // Ru
    Thermochimica::setElementMass(ctx, 46, 0.34);  // Pd

    Thermochimica::thermochimica(ctx);

    if (ctx.isSuccess()) {
        auto [nSoln, nCond] = Thermochimica::getNumberPhasesSystem(ctx);

        std::cout << "Equilibrium phases:\n\n";

        // Check each phase
        for (int i = 0; i < nSoln; ++i) {
            std::string phaseName = Thermochimica::getPhaseNameAtIndex(ctx, i);
            auto [moles, info] = Thermochimica::getMolesPhase(ctx, phaseName);

            if (info == 0 && moles > 1e-10) {
                std::cout << phaseName << ": " << moles << " mol\n";

                // Get species in this phase
                for (int j = 0; j < Thermochimica::getNumberSpeciesSystem(ctx); ++j) {
                    std::string specName = Thermochimica::getSpeciesAtIndex(ctx, j);
                    auto [molFrac, chemPot, sinfo] =
                        Thermochimica::getOutputSolnSpecies(ctx, phaseName, specName);

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

[← Back to README](../README.md) | [Databases](databases.md) | [Error Handling →](error-handling.md)
