# Thermochimica

A computational thermodynamics library written in C++17 that determines equilibrium phases and compositions for prescribed chemical conditions (temperature, pressure, composition) using Gibbs Energy Minimization (GEM).

## Features

- **Context-based architecture**: Thread-safe design with explicit state management
- **Modern C++17**: Uses Eigen for linear algebra, structured bindings, and RAII
- **Multiple solution models**: IDMX, RKMP, SUBL, SUBLM, SUBG, SUBQ, QKTO, SUBI, SUBM
- **ChemSage compatibility**: Reads standard `.dat` thermodynamic database files
- **Cross-platform**: Builds on Linux, macOS, and Windows

## Quick Start

### Prerequisites

- C++17 compiler (GCC 7+, Clang 5+, MSVC 2017+)
- CMake 3.16+
- Eigen3 and GoogleTest are automatically fetched if not found

### Building

```bash
git clone https://github.com/ORNL-CEES/thermochimica.git
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

## Usage

```cpp
#include <thermochimica/Thermochimica.hpp>

int main() {
    Thermochimica::ThermoContext ctx;

    // Load database
    Thermochimica::setThermoFilename(ctx, "data/C-O.dat");
    Thermochimica::parseCSDataFile(ctx);

    // Set conditions
    Thermochimica::setStandardUnits(ctx);
    Thermochimica::setTemperaturePressure(ctx, 1000.0, 1.0);  // 1000 K, 1 atm
    Thermochimica::setElementMass(ctx, 6, 1.0);   // 1 mol Carbon
    Thermochimica::setElementMass(ctx, 8, 2.0);   // 2 mol Oxygen

    // Calculate equilibrium
    Thermochimica::thermochimica(ctx);

    // Check results
    if (ctx.infoThermo() == 0) {
        Thermochimica::printResults(ctx);

        // Get specific values
        auto [moles, info] = Thermochimica::getSolnPhaseMol(ctx, "gas_ideal");
        if (info == 0) {
            std::cout << "Gas phase: " << moles << " mol\n";
        }
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
    Thermochimica::setThermoFilename(ctx, "data/C-O.dat");
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

## Migration from Fortran API

If migrating from the Fortran-backed C++ API, a compatibility layer is available:

```cpp
#include <thermochimica/ThermochimicaCompat.hpp>
using namespace Thermochimica::Compat;

// Old-style API calls work with minimal changes:
setThermoFilename("database.dat");
parseThermoFile();
setStandardUnits();
setTemperaturePressure(1000.0, 1.0);
setElementMass(26, 0.5);  // Fe
thermochimica();
```

See `docs/migration_guide.md` for detailed migration instructions.

## Documentation

- `docs/migration_guide.md` - Migration from Fortran API
- `doc/thermochimica_user_manual.pdf` - Theory and concepts

## License

Thermochimica has a [BSD 3-clause open-source license](LICENSE).
