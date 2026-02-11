# Thermochimica

A computational thermodynamics library written in C++17 that determines equilibrium phases and compositions for prescribed chemical conditions (temperature, pressure, composition) using Gibbs Energy Minimization (GEM).

## Features

- **Modern object-oriented API**: Clean C++ class-based interface with RAII and encapsulation
- **Thread-safe design**: Independent instances for parallel calculations
- **Modern C++17**: Uses Eigen for linear algebra, structured bindings, and move semantics
- **Multiple solution models**: IDMX, RKMP, SUBL, SUBLM, SUBG, SUBQ, QKTO, SUBI
- **ChemSage compatibility**: Reads standard `.dat` thermodynamic database files
- **Flexible units**: Temperature (K/C/F/R), pressure (atm/bar/Pa), mass (moles/grams)
- **Extensible**: Strategy pattern for custom solvers and models
- **Cross-platform**: Builds on Linux, macOS, and Windows

## Documentation

| Document | Description |
|----------|-------------|
| [Getting Started](docs/getting-started.md) | Quick start guide with minimal example |
| [ThermoClass API](docs/class_based_api.md) | Modern class-based API documentation |
| [API Reference](docs/api-reference.md) | Complete method reference |
| [Migration Guide](docs/MIGRATION_FROM_FREE_FUNCTIONS.md) | Migrate from old free-function API |
| [Building](docs/building.md) | Build instructions and project integration |
| [Configuration](docs/configuration.md) | Tolerances, units, and solver options |
| [Databases](docs/databases.md) | ChemSage format and available databases |
| [Examples](docs/examples.md) | Code examples and tutorials |
| [Error Handling](docs/error-handling.md) | Error codes and troubleshooting |
| [Architecture](docs/architecture.md) | Internal design for contributors |

### Quick Links

- **Run my first calculation** → [Getting Started](docs/getting-started.md)
- **Integrate into my project** → [Building - CMake Integration](docs/building.md#cmake-integration)
- **Look up a function** → [API Reference](docs/api-reference.md)
- **Customize solver tolerances** → [Configuration - Tolerances](docs/configuration.md#tolerances)
- **Understand an error code** → [Error Handling](docs/error-handling.md)
- **Run parallel calculations** → [Examples - Thread Safety](docs/examples.md#thread-safe-parallel-calculations)

## Quick Start

### Prerequisites

- C++17 compiler (GCC 7+, Clang 5+, MSVC 2017+)
- CMake 3.16+
- Eigen3 (auto-fetched if not found)

### Building

```bash
git clone https://github.com/dschwen/thermochimica.git
cd thermochimica
mkdir build && cd build
cmake ..
make -j
```

### Running Tests

```bash
cd build
ctest --output-on-failure
```

## Minimal Example

```cpp
#include <thermochimica/ThermoClass.hpp>

int main() {
    using namespace Thermochimica;

    // Create Thermochimica instance
    ThermoClass thermo;

    // Load database and set units
    thermo.loadDatabase("data/CO.dat");
    thermo.setStandardUnits();

    // Set conditions and composition
    thermo.setTemperaturePressure(1000.0, 1.0);
    thermo.setElementMass(6, 1.0);   // C
    thermo.setElementMass(8, 2.0);   // O

    // Run equilibrium calculation
    thermo.calculate();

    if (thermo.isSuccess()) {
        double G = thermo.getGibbsEnergy();
        // G ≈ -629533 J for CO2 at 1000 K
        thermo.printResults();
    }

    return 0;
}
```

## Thread Safety

Each `ThermoClass` instance is independent, enabling parallel calculations:

```cpp
#include <thermochimica/ThermoClass.hpp>
#include <thread>
#include <vector>

void calculate(double temperature) {
    Thermochimica::ThermoClass thermo;
    thermo.loadDatabase("data/CO.dat");
    thermo.setStandardUnits();
    thermo.setTemperaturePressure(temperature, 1.0);
    thermo.setElementMass(6, 1.0);
    thermo.calculate();
    // Results are local to this instance
}

int main() {
    std::vector<std::thread> threads;
    for (double T = 500; T <= 2000; T += 100) {
        threads.emplace_back(calculate, T);
    }
    for (auto& t : threads) {
        t.join();
    }
    return 0;
}
```

## License

Thermochimica has a [BSD 3-clause open-source license](LICENSE).
