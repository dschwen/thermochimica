# Thermochimica Interfaces

This directory contains the abstract base class interfaces for Thermochimica's object-oriented architecture.

## Overview

These interfaces define contracts for different components of the thermodynamic solver, enabling:
- **Extensibility**: New implementations can be added without modifying existing code
- **Testability**: Components can be mocked for unit testing
- **Modularity**: Components can be mixed and matched independently
- **Strategy Pattern**: Runtime selection of algorithms

## Interface Descriptions

### IThermodynamicModel.hpp
Abstract interface for thermodynamic mixing models (QKTO, RKMP, SUBL, SUBG, SUBQ, ideal).

**Key Methods:**
- `computeExcessGibbs()` - Compute excess Gibbs energy for a solution phase
- `getModelType()` - Return model type enum
- `canHandle()` - Check if model applies to a phase

**Implementations:**
- `IdealMixingModel` - Ideal solution (no excess)
- `QKTOModel` - Quasi-chemical Kohler-Toop
- `RKMPModel` - Redlich-Kister-Muggianu
- `SUBLModel` - Sublattice (Compound Energy Formalism)
- `SUBGModel` - Modified Quasichemical
- `SUBQModel` - Quasichemical variant

### ISolver.hpp
Abstract interface for GEM solver strategies.

**Key Methods:**
- `solve()` - Run solver to convergence
- `initialize()` - Initialize solver state
- `isConverged()` - Check convergence

**Implementations:**
- `StandardGEMSolver` - Unconstrained Gibbs energy minimization
- `ConstrainedGEMSolver` - Constrained minimization with phase fraction targets

### INewtonSolver.hpp
Abstract interface for Newton direction computation.

**Key Methods:**
- `computeDirection()` - Compute Newton search direction
- `handleSingular()` - Handle singular/near-singular matrices

**Implementations:**
- `NewtonSolver` - Full Hessian Newton method
- (Future: `ReducedNewtonSolver`, `RegularizedNewtonSolver`)

### ILineSearch.hpp
Abstract interface for line search strategies.

**Key Methods:**
- `search()` - Find step length along Newton direction
- `updateState()` - Update state with computed step

**Implementations:**
- `WolfeLineSearch` - Wolfe conditions (current)
- (Future: `BacktrackingLineSearch`, `AdaptiveLineSearch`)

### IParser.hpp
Abstract interface for database file parsers.

**Key Methods:**
- `parse()` - Parse database file into state
- `getSupportedExtensions()` - Get file extensions
- `canParse()` - Check if parser supports file

**Implementations:**
- `ChemSageParser` - ChemSage .dat format
- (Future: `TDBParser`, `JSONParser`)

### IConvergenceCriterion.hpp
Abstract interface for individual convergence checks.

**Key Methods:**
- `check()` - Check if criterion is satisfied
- `getName()` - Get criterion name for logging
- `getStatusMessage()` - Get detailed status
- `getTolerance()` - Get tolerance value

**Implementations:**
- `MassBalanceCheck` - Element mass balance
- `ChemicalPotentialCheck` - Chemical potential equilibrium
- `PhaseRuleCheck` - Phase rule compliance
- `DrivingForceCheck` - Unstable phase driving forces
- `MiscibilityGapCheck` - Miscibility gap detection
- `ConstraintCheck` - Phase fraction constraint satisfaction

## Usage Example

```cpp
// Create a custom thermodynamic model
class MyCustomModel : public IThermodynamicModel {
public:
    void computeExcessGibbs(ThermoState& state, GEMState& gemState, int phaseIndex) override {
        // Your custom excess Gibbs computation
    }

    Constants::PhaseType getModelType() const override {
        return Constants::PhaseType::CUSTOM;
    }

    bool canHandle(const ThermoState& state, int phaseIndex) const override {
        // Check if this model applies to the phase
        return true;
    }

    const char* getModelName() const override {
        return "MyCustomModel";
    }
};

// Register with factory
ModelFactory factory;
factory.registerModel(Constants::PhaseType::CUSTOM,
                     std::make_unique<MyCustomModel>());
```

## Design Principles

1. **Interface Segregation**: Each interface has a single, focused responsibility
2. **Dependency Inversion**: High-level modules depend on abstractions, not implementations
3. **Open/Closed**: Open for extension (new implementations), closed for modification
4. **Liskov Substitution**: Implementations are fully substitutable

## Implementation Status

- **Phase 1** ✅: Interface definitions (current)
- **Phase 2** ⏳: Model implementations
- **Phase 3** ⏳: Solver strategy implementations
- **Phase 4** ⏳: Main Thermochimica class
- **Phase 5** ⏳: Parser implementation
- **Phase 6** ⏳: Compatibility layer
- **Phase 7** ⏳: Documentation
