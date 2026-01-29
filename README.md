# Thermochimica

A computational thermodynamics library written in C++17 that determines equilibrium phases and compositions for prescribed chemical conditions (temperature, pressure, composition) using Gibbs Energy Minimization (GEM).

## Features

- **Context-based architecture**: Thread-safe design with explicit state management
- **Modern C++17**: Uses Eigen for linear algebra, structured bindings, and RAII
- **Multiple solution models**: IDMX, RKMP, SUBL, SUBLM, SUBG, SUBQ, QKTO, SUBI
- **ChemSage compatibility**: Reads standard `.dat` thermodynamic database files
- **Flexible units**: Temperature (K/C/F/R), pressure (atm/bar/Pa), mass (moles/grams)
- **Cross-platform**: Builds on Linux, macOS, and Windows

## Documentation

| Document | Description |
|----------|-------------|
| [Getting Started](docs/getting-started.md) | Quick start guide with minimal example |
| [API Reference](docs/api-reference.md) | Complete function reference |
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
#include <thermochimica/Thermochimica.hpp>

int main() {
    Thermochimica::ThermoContext ctx;

    Thermochimica::setStandardUnits(ctx);
    Thermochimica::setThermoFilename(ctx, "data/CO.dat");
    Thermochimica::parseCSDataFile(ctx);

    Thermochimica::setTemperaturePressure(ctx, 1000.0, 1.0);
    Thermochimica::setElementMass(ctx, 6, 1.0);   // C
    Thermochimica::setElementMass(ctx, 8, 2.0);   // O

    Thermochimica::thermochimica(ctx);

    if (ctx.isSuccess()) {
        double G = Thermochimica::getGibbsEnergy(ctx);
        // G ≈ -629533 J for CO2 at 1000 K
        Thermochimica::printResults(ctx);
    }

    return 0;
}
```

## Thread Safety

Each `ThermoContext` is independent, enabling parallel calculations:

```cpp
#include <thread>
#include <vector>

void calculate(double temperature) {
    Thermochimica::ThermoContext ctx;
    Thermochimica::setThermoFilename(ctx, "data/CO.dat");
    Thermochimica::parseCSDataFile(ctx);
    Thermochimica::setStandardUnits(ctx);
    Thermochimica::setTemperaturePressure(ctx, temperature, 1.0);
    Thermochimica::setElementMass(ctx, 6, 1.0);
    Thermochimica::thermochimica(ctx);
    // Results are local to this context
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
