#pragma once

#include "../ThermoContext.hpp"
#include <Eigen/Dense>
#include <vector>

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
    /// Construct the Hessian matrix (full system)
    static void constructHessian(ThermoContext& ctx,
                                 Eigen::MatrixXd& hessian,
                                 Eigen::VectorXd& rhs);

    /// Construct the Hessian matrix (reduced system with only active elements)
    static void constructHessianReduced(ThermoContext& ctx,
                                        Eigen::MatrixXd& hessian,
                                        Eigen::VectorXd& rhs,
                                        const std::vector<int>& activeElements);

    /// Handle singular matrix case
    static void handleSingular(ThermoContext& ctx,
                               Eigen::MatrixXd& hessian,
                               Eigen::VectorXd& rhs);

    /// Handle singular matrix case (reduced system)
    static void handleSingularReduced(Eigen::MatrixXd& hessian,
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

    /// Check if phase constraints are satisfied
    static bool checkPhaseConstraints(ThermoContext& ctx);

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

/// Constrained GEM solver utilities for phase fraction constraints
/// Used by augmented Lagrangian method for phase field modeling
class ConstrainedGEM {
public:
    /// Compute current element-level phase fractions for all constrained phases
    /// @param ctx The thermochimica context
    static void computePhaseElementFractions(ThermoContext& ctx);

    /// Update Lagrange multipliers after inner GEM convergence
    /// λ += ρ * (f_p - f_p^target)
    /// @param ctx The thermochimica context
    static void updateLagrangeMultipliers(ThermoContext& ctx);

    /// Add constraint penalty gradient contributions to Newton RHS
    /// @param ctx The thermochimica context
    /// @param rhs The RHS vector to modify
    /// @param activeElements List of active element indices
    static void addConstraintGradients(ThermoContext& ctx,
                                       Eigen::VectorXd& rhs,
                                       const std::vector<int>& activeElements);

    /// Add constraint penalty Hessian contributions
    /// @param ctx The thermochimica context
    /// @param hessian The Hessian matrix to modify
    /// @param activeElements List of active element indices
    static void addConstraintHessian(ThermoContext& ctx,
                                     Eigen::MatrixXd& hessian,
                                     const std::vector<int>& activeElements);

    /// Run constrained GEM iteration with penalty terms
    /// @param ctx The thermochimica context
    static void runConstrainedGEMIteration(ThermoContext& ctx);

    /// Set up phase assemblage directly from constraints
    /// For phase field applications: constrained phases are forced into assemblage,
    /// non-constrained phases are excluded entirely.
    /// @param ctx The thermochimica context
    /// @return true if assemblage was set up successfully
    static bool setupAssemblageFromConstraints(ThermoContext& ctx);
};

} // namespace Thermochimica
