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
| `ThermoClass` | Main API class (modern, canonical interface) |
| `ThermoContext` | Internal context object containing all state |
| `ThermoState` | Thermodynamic state: phases, species, results |
| `ThermoIO` | I/O parameters: T, P, composition, output flags |
| `GEMState` | GEM solver state: iteration, convergence |
| `ParserState` | Parser state during data file reading |

### Solver Flow

1. **Input** → `loadDatabase()`, `setTemperaturePressure()`, `setElementMass()`
2. **Parsing** → Database loaded and parsed automatically
3. **Calculation** → `calculate()` runs full equilibrium calculation
4. **Output** → `printResults()`, `getOutputChemPot()`, etc.

### Error Handling

`thermo.isSuccess()` checks if calculation succeeded. Use `thermo.getErrorMessage()` for descriptions.

## API Usage

The canonical API is the modern class-based `ThermoClass`:

```cpp
#include <thermochimica/ThermoClass.hpp>

int main() {
    using namespace Thermochimica;

    ThermoClass thermo;

    thermo.loadDatabase("data/C-O.dat");
    thermo.setStandardUnits();
    thermo.setTemperaturePressure(1000.0, 1.0);
    thermo.setElementMass(6, 1.0);  // Carbon

    thermo.calculate();

    if (thermo.isSuccess()) {
        thermo.printResults();
    }
    return 0;
}
```

## Legacy Code

For accessing internal state or advanced features:

```cpp
// Access underlying context if needed for advanced use
ThermoContext& ctx = thermo.getContext();

// Direct state access
double T = ctx.io->dTemperature;
auto& state = *ctx.thermo;
```

See `docs/MIGRATION_FROM_FREE_FUNCTIONS.md` if migrating from old code.

## Dependencies

- C++17 compiler - **use clang++, not g++**
- CMake 3.16+
- Eigen3 (auto-fetched if not found)
- GoogleTest (auto-fetched if not found)

## Coding Conventions

### Naming

- PascalCase for class names: `ThermoContext`, `GEMSolver`
- camelCase for functions: `setTemperaturePressure()`, `getOutputChemPot()`
- Member variables: no prefix (use `this->` for disambiguation if needed)

### Style

- Public API is through `ThermoClass` member functions
- Use Eigen for matrix/vector operations
- Prefer `std::optional` and `std::tuple` for multi-value returns
- Use structured bindings: `auto [value, info] = thermo.getResult(name);`
- Follow RAII principles for resource management
