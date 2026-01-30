#pragma once

#include <Eigen/Dense>
#include <vector>

namespace Thermochimica {

/// Mode of phase fraction constraint
enum class PhaseConstraintMode {
    None,           ///< Phase fraction is free (unconstrained)
    Fixed,          ///< Phase fraction fixed at target value
    Minimum,        ///< Phase fraction >= target (future extension)
    Maximum         ///< Phase fraction <= target (future extension)
};

/// Individual phase constraint data
struct PhaseConstraint {
    PhaseConstraintMode mode = PhaseConstraintMode::None;
    double targetFraction = 0.0;      ///< Target phase fraction [0, 1]
    double lagrangeMultiplier = 0.0;  ///< Lagrange multiplier for augmented Lagrangian
    double currentFraction = 0.0;     ///< Current computed element fraction
};

/// Phase constraint state for constrained GEM solver
/// Used for phase field modeling applications where phase fractions
/// are defined at the element level: f_p = (sum of element moles in phase p) / (total element moles)
struct PhaseConstraints {
    /// Constraints for solution phases (indexed by phase index)
    std::vector<PhaseConstraint> solnPhaseConstraints;

    /// Constraints for pure condensed phases (indexed by species index)
    std::vector<PhaseConstraint> condPhaseConstraints;

    /// Penalty parameter (rho) for augmented Lagrangian method
    /// This value grows during iteration; reset() restores it to initialPenaltyParameter
    double penaltyParameter = 1.0;

    /// Initial penalty parameter (user-configurable via setConstraintPenaltyParameter)
    /// reset() restores penaltyParameter to this value
    double initialPenaltyParameter = 1.0;

    /// Growth rate for penalty parameter each outer iteration
    double penaltyGrowthRate = 10.0;

    /// Tolerance for constraint satisfaction
    double constraintTolerance = 1e-4;

    /// Maximum outer iterations for augmented Lagrangian loop
    int maxOuterIterations = 20;

    /// Current outer iteration number
    int currentOuterIteration = 0;

    /// Default constructor
    PhaseConstraints() = default;

    /// Allocate arrays based on system dimensions
    /// @param nSolnPhasesSys Number of solution phases in system
    /// @param nConPhasesSys Number of pure condensed phase species in system
    void allocate(int nSolnPhasesSys, int nConPhasesSys);

    /// Reset constraint state (keeps constraint definitions, resets solver state)
    void reset();

    /// Clear all constraints
    void clear();

    /// Check if there are any active constraints
    bool hasActiveConstraints() const;

    /// Get number of active fixed constraints
    int getNumActiveConstraints() const;

    /// Check if all constraints are satisfied within tolerance
    bool areConstraintsSatisfied() const;

    /// Compute total constraint violation (sum of squared violations)
    double getTotalConstraintViolation() const;
};

} // namespace Thermochimica
