# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Thermochimica is a computational thermodynamics library written in C++17 that determines equilibrium phases and compositions for prescribed chemical conditions (temperature, pressure, composition) using Gibbs Energy Minimization (GEM).

## Build Commands

```bash
# Configure and build (recommended first-time setup)
mkdir build && cd build
cmake ..
make -j

# Run test suite
cd build && ctest --output-on-failure

# Build with debug symbols
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j

# Clean build
rm -rf build && mkdir build && cd build && cmake .. && make -j
```

## Architecture

### Core Design

The library uses a **context-based architecture** where all state is stored in a `ThermoContext` object. This design enables:
- Thread-safe parallel calculations
- Multiple simultaneous thermodynamic systems
- Clear ownership and lifetime management

### Directory Structure

- `include/thermochimica/` - Public headers
- `src/context/` - State management classes (ThermoContext, ThermoState, etc.)
- `src/parser/` - ChemSage data file parsing
- `src/setup/` - Input validation and initialization
- `src/models/` - Thermodynamic models (RKMP, SUBL, SUBG, SUBQ, QKTO, etc.)
- `src/solver/` - GEM solver implementation
- `src/postprocess/` - Result processing and output
- `src/api/` - Public API functions
- `tests/` - Unit and integration tests
- `data/` - Thermodynamic databases in ChemSage format

### Key Classes

| Class | Purpose |
|-------|---------|
| `ThermoContext` | Main context object containing all state |
| `ThermoState` | Thermodynamic state: phases, species, results |
| `ThermoIO` | I/O parameters: T, P, composition, output flags |
| `GEMState` | GEM solver state: iteration, convergence |
| `ParserState` | Parser state during data file reading |

### Solver Flow

1. **Input** → `setThermoFilename()`, `setTemperaturePressure()`, `setElementMass()`
2. **Parsing** → `parseCSDataFile()` reads ChemSage `.dat` databases
3. **Calculation** → `thermochimica()` runs full equilibrium calculation
4. **Output** → `printResults()`, `getOutputChemPot()`, etc.

### Error Handling

`ctx.infoThermo()` returns an integer code (0 = success). Use `getErrorMessage()` for descriptions.

## API Usage

```cpp
#include <thermochimica/Thermochimica.hpp>

int main() {
    Thermochimica::ThermoContext ctx;

    Thermochimica::setThermoFilename(ctx, "data/C-O.dat");
    Thermochimica::parseCSDataFile(ctx);

    Thermochimica::setStandardUnits(ctx);
    Thermochimica::setTemperaturePressure(ctx, 1000.0, 1.0);
    Thermochimica::setElementMass(ctx, 6, 1.0);  // Carbon

    Thermochimica::thermochimica(ctx);

    if (ctx.infoThermo() == 0) {
        Thermochimica::printResults(ctx);
    }
    return 0;
}
```

## Compatibility Layer

For migrating from the old Fortran API, use the compatibility header:

```cpp
#include <thermochimica/ThermochimicaCompat.hpp>
using namespace Thermochimica::Compat;

// Old-style API calls now work:
setThermoFilename("database.dat");
parseThermoFile();
thermochimica();
```

See `docs/migration_guide.md` for detailed migration instructions.

## Dependencies

- C++17 compiler (GCC 7+, Clang 5+, MSVC 2017+)
- CMake 3.16+
- Eigen3 (auto-fetched if not found)
- GoogleTest (auto-fetched if not found)

## Coding Conventions

### Naming

- PascalCase for class names: `ThermoContext`, `GEMSolver`
- camelCase for functions: `setTemperaturePressure()`, `getOutputChemPot()`
- Member variables: no prefix (use `this->` for disambiguation if needed)

### Style

- All public API functions take `ThermoContext&` as first parameter
- Use Eigen for matrix/vector operations
- Prefer `std::optional` and `std::tuple` for multi-value returns
- Use structured bindings: `auto [value, info] = getResult(ctx, name);`
