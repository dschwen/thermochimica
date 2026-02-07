# Bug Fixes: Pure Condensed Phase Handling

## Summary

Fixed two critical bugs in the handling of pure condensed phases (e.g., graphite) that were causing mass conservation violations and incorrect phase assemblages.

## Bugs Found and Fixed

### Bug #1: `dMolesSpecies` Not Set When Adding Pure Condensed Phases

**File**: `src/solver/PhaseAssemblage.cpp`
**Function**: `addPureConPhase()`
**Line**: 639

**Problem**: When a pure condensed phase was added to the assemblage, only `dMolesPhase` was set, but `dMolesSpecies` was never updated. This caused:
- Mass balance calculations to fail (species moles not counted)
- Inconsistent state between phase moles and species moles
- Apparent mass conservation violations

**Root Cause**: Asymmetric handling - `removePureConPhase()` correctly zeros `dMolesSpecies`, but `addPureConPhase()` didn't set it.

**Fix**: Added one line to synchronize `dMolesSpecies` with `dMolesPhase`:

```cpp
// Add to assemblage
thermo.iAssemblage(thermo.nConPhases) = speciesIndex + 1;
thermo.dMolesPhase(thermo.nConPhases) = thermo.tolerances[kTolPhaseMoles] * 10.0;
thermo.dMolesSpecies(speciesIndex) = thermo.tolerances[kTolPhaseMoles] * 10.0;  // ← ADDED
thermo.nConPhases++;
```

**Impact**:
- Restores mass conservation for pure phases
- Ensures consistent thermodynamic state
- Critical for systems with graphite, diamond, and other pure condensed species

---

### Bug #2: `getMolesPhase()` Only Searched Solution Phases

**File**: `src/api/Thermochimica.cpp`
**Function**: `getMolesPhase()`
**Lines**: 628-645

**Problem**: The `getMolesPhase()` API function only searched for **solution phases** using `getPhaseIndex()`. It never checked for **pure condensed phases**, always returning 0 for them.

**Root Cause**: The function used `getPhaseIndex()` which only finds solution phases. It should also call `getSpeciesIndex()` to find pure phases.

**Fix**: Modified the function to try both solution phases AND pure condensed phases:

```cpp
std::pair<double, int> getMolesPhase(const ThermoContext& ctx,
                                     const std::string& phaseName) {
    auto& thermo = *ctx.thermo;

    // Try as solution phase first
    int phaseIdx = thermo.getPhaseIndex(phaseName);
    if (phaseIdx >= 0) {
        // Find solution phase in assemblage
        for (int i = 0; i < thermo.nElements; ++i) {
            if (thermo.iAssemblage(i) == -(phaseIdx + 1)) {
                return {thermo.dMolesPhase(i), 0};
            }
        }
    }

    // If not found as solution phase, try as pure condensed phase
    int speciesIdx = thermo.getSpeciesIndex(phaseName);
    if (speciesIdx >= 0) {
        for (int i = 0; i < thermo.nConPhases; ++i) {
            if (thermo.iAssemblage(i) == speciesIdx + 1) {
                return {thermo.dMolesPhase(i), 0};
            }
        }
    }

    return {0.0, -1};
}
```

**Impact**:
- `getMolesPhase("C_Graphite(s)")` now works correctly
- Users can query both solution and pure phases with the same function
- Consistent API behavior

---

## Diagnostic Process

### Phase 1: Initial Discovery
- Validation tests showed discrepancies between C++ and Fortran results
- C++ missing graphite phase in Test 001 (C=1, O=1, 900K)
- Fortran correctly predicted graphite + gas assemblage

### Phase 2: Isolated Testing
Created diagnostic tool `diagnose_graphite.cpp` with 5 tests:

1. **Test A**: Verify graphite exists in database ✓
2. **Test B**: Check chemical potentials (μ_C = +12.13 J/mol → graphite should be stable) ✓
3. **Test C**: Force graphite constraint (failed - constraint not respected) ✗
4. **Test D**: Pure carbon system (1 mol C, 0 mol O) ✗ **CRITICAL**
5. **Test E**: Inspect internal state (graphite present but moles = 1e-47 ≈ 0) ✗

**Critical Finding**: Test D showed **mass conservation violation**
- Input: 1.0 mol C
- Output: 0.0 mol C (both gas and graphite had 0 moles)
- This is thermodynamically impossible!

### Phase 3: Root Cause Analysis
1. Traced through `GEMSolver::runInnerGEMLoop()`
2. Found `PhaseAssemblage::check()` was being called correctly
3. Discovered `addPureConPhase()` was missing `dMolesSpecies` update
4. Fixed bug #1

### Phase 4: API Bug Discovery
1. After fix #1, Test D still showed "Graphite: 0 mol" from `getMolesPhase()`
2. Created `diagnose_with_trace.cpp` that accessed context directly
3. Direct access showed "C_Graphite(s) = 1 mol" ✓
4. Identified that `getMolesPhase()` API was not finding pure phases
5. Fixed bug #2

## Verification

### Before Fix:
```
Test D: Pure Carbon System
  Input: C = 1.0 mol
  Output: C = 0.0 mol  ❌ MASS CONSERVATION VIOLATED
  Graphite: 0 mol
```

### After Fix:
```
Test D: Pure Carbon System
  Input: C = 1.0 mol
  Output: C = 1.0 mol  ✓ Mass conserved
  Graphite: 1 mol      ✓ Correct phase
```

### Test Suite Results:
- **All 100 integration tests pass** ✓
- **No regressions introduced** ✓
- **Validation tests now produce correct phase assemblages** ✓

## Files Modified

1. `src/solver/PhaseAssemblage.cpp` - Added `dMolesSpecies` synchronization
2. `src/api/Thermochimica.cpp` - Enhanced `getMolesPhase()` to check pure phases

## Tools Created for Diagnosis

1. `tools/diagnose_graphite.cpp` - 5 comprehensive diagnostic tests
2. `tools/diagnose_with_trace.cpp` - Detailed state inspection with mass balance check
3. `tools/check_elements.cpp` - Verify element ordering and species names
4. `tools/batch_validate.cpp` - Enhanced to include phase assemblage output
5. `tools/compare_results.py` - Compare C++ vs Fortran JSON outputs

## Impact on Validation

The fixes ensure that the C++ implementation now correctly handles:
- Pure carbon systems (graphite/diamond precipitation)
- Mixed gas + condensed phase equilibria
- Systems where pure phases should be added during iteration
- API queries for both solution and pure condensed phases

## Lessons Learned

1. **Symmetry is Critical**: If `remove` updates state variables, `add` must also update them
2. **API Completeness**: Generic query functions (`getMolesPhase`) should handle all phase types
3. **Mass Conservation Tests**: Pure component systems are excellent diagnostic cases
4. **Direct State Access**: Sometimes bypassing the API reveals bugs IN the API

## References

- Diagnostic plan: `tools/PHASE1_DIAGNOSTICS.md`
- Diagnostic results: `tools/DIAGNOSTIC_RESULTS.md`
- Validation workflow: `tools/VALIDATION_WORKFLOW.md`
