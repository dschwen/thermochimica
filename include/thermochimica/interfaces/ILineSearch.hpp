/// @file ILineSearch.hpp
/// @brief Interface for line search strategies
/// @details Defines the contract for line search along Newton direction

#pragma once

#include <Eigen/Dense>

namespace Thermochimica {

// Forward declarations
class ThermoState;
class GEMState;

/// @brief Abstract interface for line search strategies
/// @details Implementations provide different line search methods:
/// - BacktrackingLineSearch: Simple backtracking with Armijo condition
/// - WolfeLineSearch: Wolfe conditions (curvature + sufficient decrease)
/// - AdaptiveLineSearch: Adaptive step size based on convergence history
class ILineSearch {
public:
    virtual ~ILineSearch() = default;

    /// @brief Perform line search along Newton direction
    /// @param state Thermodynamic state
    /// @param gemState GEM solver state
    /// @param direction Newton direction vector
    /// @param stepLength Output: computed step length (0 < alpha <= 1)
    /// @details Finds step length that satisfies line search conditions.
    /// May modify stepLength in place or use as initial guess.
    virtual void search(ThermoState& state,
                       GEMState& gemState,
                       const Eigen::VectorXd& direction,
                       double& stepLength) = 0;

    /// @brief Update state variables with computed step
    /// @param state Thermodynamic state to update
    /// @param gemState GEM solver state to update
    /// @param direction Newton direction vector
    /// @param stepLength Step length computed by search()
    /// @details Updates element potentials, phase moles, species compositions, etc.
    virtual void updateState(ThermoState& state,
                            GEMState& gemState,
                            const Eigen::VectorXd& direction,
                            double stepLength) = 0;

    /// @brief Get line search name for logging
    /// @return Line search strategy name string
    virtual const char* getLineSearchName() const = 0;
};

} // namespace Thermochimica
