# Gibbs Energy Minimization (GEM) Solver Design Document

## Overview

The Gibbs Energy Minimization (GEM) solver computes thermodynamic equilibrium by minimizing the total Gibbs energy of a multi-component, multi-phase system subject to mass balance constraints. This document describes the algorithm and implementation in the Thermochimica C++ codebase.

### Problem Statement

Given:
- Temperature T and pressure P
- Total moles of each element (b_j for element j)
- Thermodynamic database with species and their Gibbs energy functions

Find:
- The stable phase assemblage
- Moles of each phase
- Composition (mole fractions) of each solution phase
- Element potentials (Lagrange multipliers)

Such that the total Gibbs energy G is minimized:

```
G = Σᵢ nᵢ μᵢ   (sum over all species i)
```

Subject to mass balance constraints:

```
Σᵢ nᵢ aᵢⱼ = bⱼ   for each element j
```

Where:
- nᵢ = moles of species i
- μᵢ = chemical potential of species i
- aᵢⱼ = stoichiometric coefficient of element j in species i
- bⱼ = total moles of element j

---

## Algorithm Architecture

The solver consists of four main stages executed in sequence:

```
┌─────────────────────────────────────────────────────────┐
│                    thermochimica()                       │
├─────────────────────────────────────────────────────────┤
│  1. init()           - Initialize context               │
│  2. checkSystem()    - Validate inputs                  │
│  3. compThermoData() - Compute Gibbs energies at T      │
│  4. setup()          - Initial phase assemblage         │
│     ├── levelingSolver()     - Linear estimation        │
│     └── postLevelingSolver() - Mole fraction estimates  │
│  5. solve()          - GEM iterations                   │
│     └── GEMSolver::solve()                              │
│  6. postProcess()    - Compute output quantities        │
└─────────────────────────────────────────────────────────┘
```

---

## Stage 1: Leveling Solver (Initial Estimation)

**File:** `src/solver/LevelingSolver.cpp`

### Purpose

The leveling solver provides an initial estimate of the phase assemblage and element potentials. It uses the technique of Eriksson and Thompson (CALPHAD V.13, 1989) which temporarily treats all species as pure stoichiometric phases, making the problem linear.

### Algorithm

```
1. SELECT initial assemblage:
   - For each element, select the pure condensed species with lowest
     Gibbs energy that contains that element

2. ITERATE until converged (max 100 iterations):
   a. Solve for element potentials λⱼ:
      - Build stoichiometry matrix A from current assemblage
      - Solve: A^T λ = G°/RT  (normalized Gibbs energies)

   b. Compute driving force for ALL species:
      - DF_i = G°ᵢ/RT - Σⱼ λⱼ aᵢⱼ

   c. CHECK convergence:
      - If min(DF_i) ≥ -tolerance → CONVERGED

   d. SWAP phases:
      - Find species with most negative driving force
      - Replace least favorable species in assemblage

3. COMPUTE phase moles:
   - Solve: A · n = b  (mass balance)
```

### Key Data Structures

```cpp
iAssemblage[e]      // Species index in position e (1-based, or negative for solution phases)
dElementPotential[j] // Lagrange multiplier for element j (λⱼ)
dMolesPhase[e]      // Moles of phase in position e
```

### Driving Force Interpretation

- **Negative driving force**: Species is more stable than current potentials allow → should be added
- **Positive driving force**: Species is less stable → should be removed or not added
- **Zero driving force**: Species is in equilibrium with current potentials

---

## Stage 2: Post-Leveling (Solution Phase Initialization)

**File:** `src/solver/LevelingSolver.cpp` (function `postLevelingSolver`)

### Purpose

Estimates initial mole fractions for solution phase species based on element potentials from the leveling solver.

### Algorithm

For each solution phase:
```
1. For each species i in the phase:
   - Compute μ*ᵢ = Σⱼ λⱼ aᵢⱼ / particles_per_mole
   - Compute xᵢ ∝ exp(μ*ᵢ - G°ᵢ)

2. Normalize: xᵢ = xᵢ / Σₖ xₖ
```

---

## Stage 3: GEM Solver (Main Iteration Loop)

**File:** `src/solver/GEMSolver.cpp`

### Main Loop Structure

```cpp
for (iter = 1; iter <= MAX_ITER; ++iter) {
    // 1. Compute Newton direction
    if (nSolnPhases > 0 || functionNorm > tolerance) {
        GEMNewton::compute(ctx);
        GEMLineSearch::search(ctx);
    }

    // 2. Check and adjust phase assemblage
    PhaseAssemblage::check(ctx);

    // 3. Check convergence
    if (ConvergenceChecker::check(ctx)) {
        break;  // SUCCESS
    }
}
```

### Pure Condensed Phase Shortcut

For systems with only pure condensed phases (no solution phases), the problem is linear. If the initial assemblage from leveling satisfies mass balance within tolerance, the solver exits immediately.

---

## Stage 4: Newton Direction Computation

**File:** `src/solver/GEMNewton.cpp`

### Purpose

Computes the search direction for the next iteration by solving a linearized system of equations (KKT conditions for the constrained optimization).

### System Structure

The Newton system has the form:

```
[ H   A^T ] [ Δλ ]   [ r_λ ]
[ A   0   ] [ Δn ] = [ r_n ]
```

Where:
- H = Hessian contribution from solution phases
- A = Stoichiometry matrix
- Δλ = Update to element potentials
- Δn = Update to phase moles
- r = Residual vector

### Hessian Construction

The Hessian matrix has three main blocks:

**Block 1: Element-Element (from solution phases)**
```cpp
H[i,j] += Σₖ (nₖ · aₖᵢ · aₖⱼ) / particles²
```
Where the sum is over species k in solution phases.

**Block 2: Element-Solution Phase**
```cpp
H[j, nElements + p] = Σₖ (xₖ · aₖⱼ) / particles
```
Where the sum is over species k in solution phase p.

**Block 3: Element-Condensed Phase**
```cpp
H[j, nElements + nSolnPhases + p] = aᵢⱼ / particles
```
Where i is the species index for pure condensed phase p.

### Right-Hand Side (Residual)

**Element constraints:**
```cpp
r[j] = b[j] - Σᵢ (nᵢ · aᵢⱼ) + correction_terms
```

**Solution phase constraints:**
```cpp
r[nElements + p] = G_excess_phase[p]
```

**Pure condensed phase constraints:**
```cpp
r[nElements + nSolnPhases + p] = G°[species_p]
```

### Linear Solve

Uses Eigen's `PartialPivLU` (equivalent to LAPACK DGESV):
```cpp
Eigen::PartialPivLU<Eigen::MatrixXd> lu(hessian);
dUpdateVar = lu.solve(rhs);
```

---

## Stage 5: Line Search

**File:** `src/solver/GEMLineSearch.cpp`

### Purpose

Finds an appropriate step length along the Newton direction to ensure:
1. Sufficient decrease in the objective (mass balance residual)
2. Phase moles remain non-negative
3. Convergence is not disrupted

### Algorithm

```
1. INITIALIZE step length γ = 1.0

2. CONSTRAIN step based on:
   - Maximum element potential change (prevents large jumps)
   - Phase moles going negative (prevents infeasible steps)
   - Iteration count (more conservative at high iterations)

3. WOLFE LINE SEARCH (max 5 iterations):
   a. Save current state
   b. Apply update: x_new = x_old + γ · Δx
   c. Compute new function norm (mass balance residual)
   d. If sufficient decrease → ACCEPT
   e. Else: γ = γ / 2, REPEAT

4. UPDATE tracked quantities for next iteration
```

### Step Length Constraints

```cpp
// Limit by element potential change
for each element j:
    γ = min(γ, maxGamma / |Δλⱼ|)

// Prevent negative phase moles
for each phase p:
    if (n_new < 0 and n_current > 0):
        γ = min(γ, 0.9 · n_current / |n_new - n_current|)
```

### Wolfe Conditions (Simplified)

Accept step if:
- Function norm < 10⁻⁶ (excellent progress), OR
- Relative decrease > 0.1% (sufficient decrease), OR
- Plateau detected (relative change ~1.0)

---

## Stage 6: Phase Assemblage Management

**File:** `src/solver/PhaseAssemblage.cpp`

### Purpose

Dynamically adjusts which phases are included in the equilibrium calculation based on their stability (driving force).

### Phase Removal

Phases are removed when:
- Phase moles become negative or below tolerance
- Driving force becomes positive (phase is unstable)

```cpp
if (dMolesPhase[p] < tolerance) {
    removePhase(p);
}
```

### Phase Addition

Phases are added when:
- Driving force is significantly negative (phase would lower total Gibbs energy)
- Gibbs phase rule allows (P ≤ C components)

```cpp
// Compute driving force for unstable solution phase
DF = Σᵢ xᵢ · (μᵢ - μ*ᵢ)

if (DF > tolerance) {
    addSolnPhase(phase);
}
```

### Phase Swapping

When the phase rule prevents adding a new phase, the solver may swap:
1. Find the pure condensed phase with worst (most positive) driving force
2. If new phase has better driving force, remove the pure phase and add the solution phase

### Assemblage Indexing Convention

```
iAssemblage[i] > 0  : Pure condensed phase, value = species index (1-based)
iAssemblage[i] < 0  : Solution phase, value = -(phaseIndex + 1)
iAssemblage[i] = 0  : Empty slot
```

---

## Stage 7: Convergence Checking

**File:** `src/solver/CheckConvergence.cpp`

### Convergence Criteria

The solver checks multiple criteria in order:

1. **Function Norm**: Mass balance residual < tolerance
2. **Phase Validity**: All active slots have assigned phases
3. **Phase Moles**: All phase moles ≥ 0
4. **Phase Rule**: P ≤ C (number of phases ≤ number of components)
5. **Mass Balance**: |Σᵢ nᵢ aᵢⱼ - bⱼ| / max(1, bⱼ) < tolerance for all j
6. **Chemical Potential**: Element potentials are finite (no NaN/Inf)
7. **Site Fractions**: Sum to 1.0 for each sublattice (sublattice phases only)
8. **Unstable Phases**: No unstable phase has favorable driving force
9. **Miscibility**: No undetected miscibility gaps (placeholder)

### Tolerance Values

```cpp
kTolFunctionNorm   = 1e-8   // Mass balance residual norm
kTolMassBalance    = 1e-6   // Relative mass balance error
kTolDrivingForce   = 1e-6   // Minimum driving force to add phase
kTolPhaseMoles     = 1e-15  // Minimum moles for phase to exist
kTolSiteFraction   = 1e-6   // Site fraction sum deviation from 1
```

---

## Data Flow Diagram

```
┌──────────────┐     ┌──────────────┐     ┌──────────────┐
│   Parser     │────▶│ ThermoState  │◀───▶│  GEMState    │
│ (database)   │     │ (species,    │     │ (iteration   │
│              │     │  phases)     │     │  state)      │
└──────────────┘     └──────────────┘     └──────────────┘
                            │
                            ▼
                     ┌──────────────┐
                     │   IOState    │
                     │ (T, P, mass, │
                     │  results)    │
                     └──────────────┘
```

### Key State Variables

**ThermoState** (persistent across iterations):
- `nSpecies`, `nElements` - System dimensions
- `dStoichSpecies[i,j]` - Stoichiometry matrix
- `dStdGibbsEnergy[i]` - Standard Gibbs energies at T
- `dMolesSpecies[i]` - Moles of each species
- `dMolFraction[i]` - Mole fractions in solution phases
- `dElementPotential[j]` - Lagrange multipliers
- `iAssemblage[e]` - Current phase assemblage

**GEMState** (solver working state):
- `iterGlobal` - Current iteration number
- `dGEMFunctionNorm` - Mass balance residual norm
- `dUpdateVar` - Newton direction vector
- `lConverged` - Convergence flag
- `lSolnPhases[p]` - Solution phase stability flags
- `dDrivingForceSoln[p]` - Driving forces for unstable phases

---

## Performance Considerations

### Complexity

- **Leveling Solver**: O(N³) per iteration for N elements (matrix solve)
- **GEM Newton**: O(M³) per iteration for M = N_elements + N_phases
- **Line Search**: O(N_species × N_elements) for function evaluation

### Convergence Behavior

- **Pure phases only**: Often converges in 1-10 iterations
- **With solution phases**: Typically 50-500 iterations
- **Complex systems**: May require 1000+ iterations

### Memory Layout

Arrays use Eigen's column-major storage for efficient linear algebra:
```cpp
Eigen::MatrixXd hessian;  // Column-major by default
Eigen::VectorXd rhs;
```

---

## References

1. Eriksson, G. and Thompson, W.T. "A Procedure to Estimate Equilibrium Concentrations in Multicomponent Systems and Related Applications." CALPHAD V.13, N.4, pp.389-400 (1989).

2. Piro, M.H.A. "Computation of Thermodynamic Equilibria Pertinent to Nuclear Materials in Multi-Physics Codes." PhD Dissertation, Royal Military College of Canada (2011).

3. Hillert, M. "The Compound Energy Formalism." Journal of Alloys and Compounds 320, pp.161-176 (2001).

---

## File Index

| File | Purpose |
|------|---------|
| `GEMSolver.cpp` | Main solver loop and initialization |
| `GEMNewton.cpp` | Newton direction computation |
| `GEMLineSearch.cpp` | Line search with Wolfe conditions |
| `LevelingSolver.cpp` | Initial phase assemblage estimation |
| `PhaseAssemblage.cpp` | Dynamic phase addition/removal |
| `CheckConvergence.cpp` | Convergence criteria checking |
| `GEMSolver.hpp` | Class declarations and interfaces |
