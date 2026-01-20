# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Thermochimica is a computational thermodynamics library written in Fortran 90+ that determines equilibrium phases and compositions for prescribed chemical conditions (temperature, pressure, composition) using Gibbs Energy Minimization (GEM).

## Build Commands

```bash
# Build library and all tests (recommended first-time setup)
make test

# Run test suite
./run_tests

# Build library only
make

# Build with debug symbols (-O0 -g)
make debug

# Clean build
make clean

# Full clean (removes obj/ and bin/)
make veryclean

# Generate Doxygen documentation
make doc
```

### Running Individual Tests

```bash
# After building, run specific test
./bin/TestThermo01

# Create and run custom test
cp test/Thermo.F90 test/demo.F90
# Edit demo.F90
make
./bin/demo
```

### Running Input Scripts

```bash
./bin/InputScriptMode inputs/demo.ti
```

## Architecture

### Core Solver Flow

1. **Input** → `Thermochimica()` entry point validates T, P, composition via `CheckThermoInput()`
2. **Parsing** → `ParseCSDataFile()` reads ChemSage `.dat` thermodynamic databases from `data/`
3. **Calculation** → `CompThermoData()` computes Gibbs energies
4. **Equilibrium** → `GEMSolver()` performs iterative optimization (Newton direction + line search)
5. **Output** → `PostProcessThermo()` and optionally `WriteJSON(.TRUE.)` for JSON output

### Key Directories

- `src/module/` - Core modules with shared state (`ModuleThermo.f90` contains 100+ global variables)
- `src/gem/` - GEM solver implementation (58 files)
- `src/parser/` - ChemSage data file parsing
- `src/api/` - Public API for coupling (30 files with getter/setter functions)
- `src/` - C/C++ interface bindings (`Thermochimica-c.C`, `Thermochimica-cxx.C`)
- `test/daily/` - Unit tests (`TestThermo[01-90].F90`)
- `data/` - Thermodynamic databases in ChemSage format

### Key Modules

| Module | Purpose |
|--------|---------|
| `ModuleThermo.f90` | Core state: phase/species data, results |
| `ModuleThermoIO.f90` | I/O: input T,P, composition, output flags |
| `ModuleGEMSolver.f90` | GEM solver state: iteration, convergence |
| `ModuleParseCS.f90` | Parser state during data file reading |

### Error Handling

`INFOThermo` integer code (0 = success). Error codes defined in `src/debug/ThermoDebug.f90`.

## Coding Conventions

### Naming Prefixes

- `d` - double precision real (e.g., `dTemperature`, `dMolesPhase`)
- `i` - integer (e.g., `iAssemblage`, `nConPhases`)
- `c` - character/string (e.g., `cElementName`, `cSpeciesName`)
- `l` - logical (e.g., `lConverged`, `lDebugMode`)

### Style

- One major subroutine per `.f90` file
- PascalCase for subroutine names (e.g., `GEMSolver`, `CheckPhaseAssemblage`)
- Doxygen-style comments (`!>`, `!!`) with header blocks documenting purpose, author, date

### Test Pattern

```fortran
program TestThermo01
    USE ModuleThermoIO
    implicit none
    dTemperature = 300D0
    dPressure = 1D0
    dElementMass(6) = 1D0  ! Carbon
    call ParseCSDataFile('data/C-O.dat')
    call Thermochimica
    if (INFOThermo == 0) then
        print *, 'TestThermo01: PASS'
        call EXIT(0)
    else
        print *, 'TestThermo01: FAIL'
        call EXIT(1)
    end if
end program
```

## Dependencies

- Fortran 90+ compiler (gfortran)
- LAPACK/BLAS (macOS uses Accelerate framework)
- g++ for C/C++ interfaces

## Platform Notes

- Linux: Links `-llapack -lblas -lgfortran`
- macOS: Uses `-framework Accelerate`
- Compiler flags: `-Wall -O2 -ffree-line-length-none -fbounds-check -ffpe-trap=zero -cpp`
