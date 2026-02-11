# Phase 3 Implementation Plan: Solver Strategies

## Overview

Phase 3 involves converting solver components from static utility classes to strategy implementations. This is the most complex phase as solver components are tightly coupled and interdependent.

## Challenges Identified

1. **Tight Coupling**: GEMSolver, GEMNewton, GEMLineSearch, and PhaseAssemblage are highly interdependent
2. **State Management**: Current design uses ThermoContext; new design needs ThermoState + GEMState separation
3. **Legacy Compatibility**: Must maintain 100% backward compatibility during transition
4. **Test Coverage**: All 79 tests must continue passing

## Recommended Approach: Incremental Refactoring

Rather than attempting a complete rewrite, use an incremental approach:

### Step 1: Create Adapter Classes (Bridge Pattern)
Create thin wrapper classes that implement interfaces but delegate to existing static methods:

```cpp
class WolfeLineSearch : public ILineSearch {
    void search(...) override {
        // Reconstruct ThermoContext from state/gemState
        // Call existing GEMLineSearch::search()
    }
};
```

### Step 2: Gradual Migration
Once adapters are working:
1. Tests verify identical behavior
2. Incrementally move implementation from static methods to instance methods
3. Each small change followed by test verification

### Step 3: Refactor Internals
After external interface is stable:
1. Refactor internal implementation details
2. Remove ThermoContext dependencies
3. Optimize for the new architecture

## Files to Create

### Headers
- `include/thermochimica/solver/WolfeLineSearch.hpp`
- `include/thermochimica/solver/NewtonSolver.hpp`
- `include/thermochimica/solver/StandardGEMSolver.hpp`
- `include/thermochimica/solver/ConstrainedGEMSolver.hpp`
- `include/thermochimica/solver/PhaseAssemblageManager.hpp`

### Implementations
- `src/solver/WolfeLineSearch.cpp`
- `src/solver/NewtonSolver.cpp`
- `src/solver/StandardGEMSolver.cpp`
- `src/solver/ConstrainedGEMSolver.cpp`
- `src/solver/PhaseAssemblageManager.cpp`

## Implementation Order

1. **WolfeLineSearch** (simplest, ~250 LOC)
2. **NewtonSolver** (medium complexity, ~450 LOC)
3. **PhaseAssemblageManager** (composition, ~730 LOC)
4. **StandardGEMSolver** (complex, ~650 LOC)
5. **ConstrainedGEMSolver** (complex, ~700 LOC)

## Risk Mitigation

1. **Preserve Legacy Code**: Keep all existing files unchanged initially
2. **Parallel Implementation**: New classes coexist with old static classes
3. **Comprehensive Testing**: Run full test suite after each component
4. **Gradual Migration**: Switch one component at a time

## Estimated Effort

- **Full Phase 3**: 2 weeks (as per original plan)
- **Minimal Phase 3** (adapters only): 2-3 days
- **Current Status**: Interfaces defined, implementation pending

## Next Steps

Given project complexity and time investment so far, recommend:

**Option A: Minimal Phase 3**
- Create adapter classes only
- Verify tests pass
- Defer full refactoring to future work
- **Benefit**: Quick progress to Phase 4 (main Thermochimica class)

**Option B: Complete Phase 3**
- Full refactoring of all solver components
- Move implementation into strategy classes
- Comprehensive internal restructuring
- **Benefit**: Clean architecture throughout
- **Cost**: 2 weeks development time

**Option C: Skip to Phase 4**
- Implement Thermochimica class using existing solver components
- Come back to Phase 3 later if needed
- **Benefit**: Working OO API sooner
- **Trade-off**: Mixed architecture (some OO, some static)

## Recommendation

**Option A (Minimal Phase 3)** provides best balance of:
- Progress toward complete OO architecture
- Risk mitigation (minimal changes)
- Quick validation (tests pass quickly)
- Foundation for future work

After Phase 4 (main class) is complete, can return to fully refactor Phase 3 internals.
