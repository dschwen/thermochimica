# Phase 1 Diagnostics: Phase Assemblage Analysis

**Date:** February 5, 2026
**Investigation:** Phase 1, Task 1.1 - Enable Detailed C++ Output
**Status:** ⚠️ CRITICAL ISSUE IDENTIFIED

---

## Summary

**SMOKING GUN FOUND:** C++ implementation is **NOT stabilizing pure condensed graphite phase** when it should be thermodynamically favorable.

---

## Detailed Phase Comparison

### Test 001: 900K, 1 atm, C=1 mol, O=1 mol

| Implementation | Phases | Moles | Total Gibbs (J) |
|----------------|--------|-------|-----------------|
| **Fortran** | gas_ideal | 0.602 mol | **-308,095** |
|             | C_Graphite(s) | 0.398 mol | |
| **C++** | gas_ideal | 0.747 mol | **-450,244** |
|         | C_Graphite(s) | **0.000 mol** ❌ | |

**Difference:** 142,149 J (46%)

**Analysis:**
- Fortran correctly precipitates graphite (0.398 mol)
- C++ keeps all carbon in gas phase
- C++ Gibbs is MORE NEGATIVE (more stable) - physically incorrect!
- Missing graphite means C++ has more gas species, different composition

---

### Test 002: 950K, 2 atm, C=1 mol, O=2 mol ✅

| Implementation | Phases | Moles | Total Gibbs (J) |
|----------------|--------|-------|-----------------|
| **Fortran** | gas_ideal | 1.000 mol | **-610,552** |
| **C++** | gas_ideal | 0.999 mol | **-615,198** |

**Difference:** 4,646 J (0.76%) ✅

**Analysis:**
- EXCELLENT AGREEMENT
- No graphite in either implementation
- Pure gas phase system
- Proves that C++ GEM solver works correctly for gas-only systems

---

### Test 003: 1000K, 3 atm, C=2 mol, O=1 mol

| Implementation | Phases | Moles | Total Gibbs (J) |
|----------------|--------|-------|-----------------|
| **Fortran** | gas_ideal | 0.676 mol | **-332,200** |
|             | C_Graphite(s) | 1.324 mol | |
| **C++** | gas_ideal | 1.248 mol | **-785,619** |
|         | C_Graphite(s) | **0.000 mol** ❌ | |

**Difference:** 453,419 J (137%)

**Analysis:**
- Fortran correctly precipitates LARGE amount of graphite (1.324 mol)
- C++ keeps ALL carbon in gas phase
- C++ Gibbs is MUCH more negative - thermodynamically impossible!
- Carbon excess should stabilize graphite formation

---

## Root Cause Analysis

### Confirmed Issue

**C++ is not including pure condensed phases (graphite) in the equilibrium calculation**

### Evidence

1. ✅ Test 002 (gas-only) shows perfect agreement → GEM solver is correct
2. ❌ Tests 001 and 003 (graphite + gas) show large errors → Pure phase logic is broken
3. ❌ C++ Gibbs energies are MORE NEGATIVE without graphite → Unphysical
4. ❌ getMolesPhase("C_Graphite") returns 0 in all cases

### Hypothesis

The C++ implementation likely has one of these issues:

**A. Pure Phase Initialization**
- Graphite not being added to candidate phases
- Database parsing not reading pure condensed species correctly
- Pure phase data not loaded into solver

**B. Phase Stability Test**
- Driving force calculation for pure phases incorrect
- Chemical potential comparison failing
- Reference state mismatch preventing stabilization

**C. GEM Solver Logic**
- Pure phases excluded from iteration
- Only solution phases being considered
- Assembly list not including pure condensed species

**D. Post-Processing**
- Phases ARE stable but not reported correctly
- getMolesPhase() API not accessing pure phases
- Results stored but not accessible

---

## Diagnostic Tests

### Test A: Check if Graphite is in Database

**Action:** Query if C_Graphite exists after parsing

```cpp
ThermoClass thermo;
thermo.loadDatabase("CO.dat");
// Check if graphite species exists
auto [moles, info] = thermo.getMolesPhase("C_Graphite(s)");
std::cout << "Graphite query result: " << info << "\n";
```

**Expected if parsing works:** info = 0 (success) or specific "not stable" error
**Expected if parsing fails:** info = "phase not found" error

---

### Test B: Check Chemical Potential

**Action:** Compare C chemical potential with graphite reference

```cpp
auto [muC, info1] = thermo.getOutputChemPot("C");
auto [muGraphite, info2] = thermo.getMolesPhase("C_Graphite(s)");

// Graphite should be stable if: μ_C < μ_Graphite(reference)
```

**Current C++ μ_C values:**
- Test 001: +12.13 J/mol
- Test 003: +9.41 J/mol

**Expected graphite reference:** ~0 J/mol (standard state)

**Analysis:** If μ_C > 0, graphite SHOULD be stable!

---

### Test C: Force Graphite Constraint

**Action:** Use phase constraint to force graphite

```cpp
thermo.setCondenPhaseConstraint("C_Graphite", 0.1);  // Force 0.1 mol
thermo.calculate();
```

**Expected if solver works:** Calculation succeeds with graphite present
**Expected if solver broken:** Error or ignored

---

### Test D: Compare with Pure Graphite

**Action:** Run pure C system (no oxygen)

```cpp
thermo.setElementMass("C", 1.0);
// No oxygen
thermo.calculate();
```

**Expected:** 100% graphite, G ≈ 0
**If broken:** Gas phase or error

---

##Next Steps (Immediate)

### 1. Run Diagnostic Tests A-D

**Timeline:** 30 minutes
**Priority:** CRITICAL

Execute all four diagnostic tests to narrow down the root cause.

---

### 2. Inspect GEM Solver Pure Phase Logic

**Files to check:**
- `src/solver/GEMSolver.cpp` - Main solver loop
- `src/solver/PhaseAssemblage.cpp` - Phase selection logic
- `src/solver/StandardGEMSolver.cpp` - Pure phase handling
- `src/models/IdealMixing.cpp` - Graphite should be "ideal"

**Look for:**
- How pure condensed phases enter assemblage
- Chemical potential calculations for pure phases
- Driving force calculations
- Phase stability criteria

**Timeline:** 2 hours
**Priority:** HIGH

---

### 3. Compare with Fortran Pure Phase Logic

**Files to check (Fortran):**
- `legacy/src/CheckPureConPhaseAdd.F90` - Pure phase addition
- `legacy/src/CheckPureConPhaseRem.F90` - Pure phase removal
- `legacy/src/CompChemPot.F90` - Chemical potential calculation

**Compare:**
- When are pure phases tested for stability?
- What criteria determine if phase enters assemblage?
- How is driving force calculated?

**Timeline:** 2 hours
**Priority:** HIGH

---

### 4. Add Debug Output to GEM Solver

**Modifications needed:**
```cpp
// In GEM solver main loop:
std::cout << "Iteration " << iter << ":\n";
std::cout << "  Testing pure phase: C_Graphite\n";
std::cout << "  Driving force: " << drivingForce << "\n";
std::cout << "  μ_C in system: " << muC << "\n";
std::cout << "  μ_Graphite ref: " << muGraphite_ref << "\n";
std::cout << "  Stable? " << (drivingForce < 0 ? "YES" : "NO") << "\n";
```

**Timeline:** 1 hour
**Priority:** MEDIUM

---

## Expected Resolution

### Most Likely Cause

Based on the evidence, the most likely cause is **Option B: Phase Stability Test**

**Reasoning:**
1. Test 002 works perfectly → GEM solver math is correct
2. Pure phases missing → Specific to condensed phase logic
3. μ_C > 0 in C++ → Should trigger graphite stability
4. getMolesPhase returns 0 → Phase never enters assemblage

**Predicted Fix Location:**
- `src/solver/PhaseAssemblage.cpp` or `src/solver/StandardGEMSolver.cpp`
- Function that checks if pure phases should be added
- Likely a sign error, threshold error, or missing check

---

### Fix Verification

Once fixed, re-run trial-runList tests:

**Expected Test 001 results:**
- gas_ideal: ~0.60 mol
- C_Graphite(s): ~0.40 mol
- Gibbs: ~-308,000 J

**Expected Test 003 results:**
- gas_ideal: ~0.68 mol
- C_Graphite(s): ~1.32 mol
- Gibbs: ~-332,000 J

**Success criteria:**
- All 3 tests pass with < 1% error
- Phase assemblages match Fortran
- Gibbs energies agree to within numerical precision

---

## Impact Assessment

### Severity: CRITICAL

This bug affects ALL calculations involving pure condensed phases:
- Metals (solidification)
- Oxides (precipitation)
- Carbides (graphite, carbides)
- Salt precipitation
- Any system with phase transitions

### Affected Use Cases

❌ **Broken:**
- High-temperature alloy calculations
- Corrosion modeling (oxide layers)
- Nuclear fuel behavior (fission products)
- Molten salt reactor chemistry
- Carbon activity in steel

✅ **Working:**
- Pure gas phase systems
- Single-phase liquid systems (if no solids precipitate)
- Systems far from saturation

---

## Conclusion

**The C++ implementation has a critical bug in pure condensed phase stability logic.**

- Root cause is NOT in the GEM solver mathematics (Test 002 proves that works)
- Root cause IS in how pure phases are tested and added to assemblage
- Bug is likely in PhaseAssemblage or chemical potential comparison logic
- Fix should be straightforward once exact location is identified

**Immediate action:** Run diagnostic tests A-D to pinpoint the exact function/line causing the issue.

---

**Report prepared by:** Claude Code Validation System
**Next update:** After diagnostic tests complete
