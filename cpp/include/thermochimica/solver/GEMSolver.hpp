#pragma once

#include "../ThermoContext.hpp"

namespace Thermochimica {

/// Gibbs Energy Minimization solver
/// Main solver class converted from GEMSolver.f90
class GEMSolver {
public:
    /// Run the GEM solver
    /// @param ctx The thermochimica context
    /// @return Error code (0 = success)
    static int solve(ThermoContext& ctx);

    /// Initialize the GEM solver
    static void init(ThermoContext& ctx);

    /// Check if system with only pure condensed phases is converged
    static void checkSysOnlyPureConPhases(ThermoContext& ctx);

private:
    /// Single iteration of the GEM solver
    static void iterate(ThermoContext& ctx);
};

/// Newton direction computation
/// Converted from GEMNewton.f90
class GEMNewton {
public:
    /// Compute the Newton direction vector
    /// @param ctx The thermochimica context
    /// @return LAPACK INFO code (0 = success)
    static int compute(ThermoContext& ctx);

private:
    /// Construct the Hessian matrix
    static void constructHessian(ThermoContext& ctx,
                                 Eigen::MatrixXd& hessian,
                                 Eigen::VectorXd& rhs);

    /// Handle singular matrix case
    static void handleSingular(ThermoContext& ctx,
                               Eigen::MatrixXd& hessian,
                               Eigen::VectorXd& rhs);
};

/// Line search along Newton direction
/// Converted from GEMLineSearch.f90
class GEMLineSearch {
public:
    /// Perform line search
    /// @param ctx The thermochimica context
    static void search(ThermoContext& ctx);

private:
    /// Initialize line search parameters
    static void init(ThermoContext& ctx, double& stepLength);

    /// Apply Wolfe conditions
    static bool checkWolfeConditions(ThermoContext& ctx, double functionNorm);

    /// Update element potentials
    static void updateElementPotentials(ThermoContext& ctx, double stepLength);

    /// Update phase moles
    static void updatePhaseMoles(ThermoContext& ctx, double stepLength);

    /// Update species moles and mole fractions
    static void updateSpecies(ThermoContext& ctx);
};

/// Phase assemblage management
/// Converted from CheckPhaseAssemblage.f90
class PhaseAssemblage {
public:
    /// Check and adjust phase assemblage
    static void check(ThermoContext& ctx);

    /// Add a solution phase
    static bool addSolnPhase(ThermoContext& ctx, int phaseIndex);

    /// Remove a solution phase
    static bool removeSolnPhase(ThermoContext& ctx, int phaseIndex);

    /// Add a pure condensed phase
    static bool addPureConPhase(ThermoContext& ctx, int speciesIndex);

    /// Remove a pure condensed phase
    static bool removePureConPhase(ThermoContext& ctx, int speciesIndex);

    /// Swap phases
    static void swapPhases(ThermoContext& ctx);

    /// Revert to previous assemblage
    static void revert(ThermoContext& ctx);
};

/// Leveling solver for initial phase assemblage
/// Converted from LevelingSolver.f90
void levelingSolver(ThermoContext& ctx);

/// Post-leveling adjustments
void postLevelingSolver(ThermoContext& ctx);

/// Convergence checking
/// Converted from CheckConvergence.f90
class ConvergenceChecker {
public:
    /// Check all convergence criteria
    /// @param ctx The thermochimica context
    /// @return true if converged
    static bool check(ThermoContext& ctx);

private:
    /// Check mass balance residuals
    static bool checkMassBalance(ThermoContext& ctx);

    /// Check chemical potential residuals
    static bool checkChemicalPotential(ThermoContext& ctx);

    /// Check phase rule
    static bool checkPhaseRule(ThermoContext& ctx);

    /// Check site fractions (for sublattice phases)
    static bool checkSiteFractions(ThermoContext& ctx);

    /// Check unstable phases
    static bool checkUnstablePhases(ThermoContext& ctx);

    /// Check miscibility gaps
    static bool checkMiscibility(ThermoContext& ctx);
};

} // namespace Thermochimica
