/// @file INewtonSolver.hpp
/// @brief Interface for Newton direction computation
/// @details Defines the contract for computing Newton search directions

#pragma once

#include <Eigen/Dense>

namespace Thermochimica {

// Forward declarations
class ThermoState;
class GEMState;

/// @brief Abstract interface for Newton direction computation
/// @details Implementations handle Hessian construction, RHS vector assembly,
/// and solving the linear system for the Newton direction. Different
/// implementations may use full Hessian, reduced Hessian, or regularization.
class INewtonSolver {
public:
    virtual ~INewtonSolver() = default;

    /// @brief Compute Newton direction vector
    /// @param state Thermodynamic state
    /// @param gemState GEM solver state with Hessian and RHS
    /// @param direction Output: Newton direction vector
    /// @return Error code (0 = success, >0 = singular/error)
    virtual int computeDirection(ThermoState& state,
                                GEMState& gemState,
                                Eigen::VectorXd& direction) = 0;

    /// @brief Handle singular matrix cases
    /// @param hessian Hessian matrix (may be modified for regularization)
    /// @param rhs RHS vector (may be modified)
    /// @details Called when matrix is singular or near-singular.
    /// Implementations may add regularization, remove variables, etc.
    virtual void handleSingular(Eigen::MatrixXd& hessian,
                               Eigen::VectorXd& rhs) = 0;

    /// @brief Get solver name for logging
    /// @return Newton solver name string
    virtual const char* getNewtonSolverName() const = 0;
};

} // namespace Thermochimica
