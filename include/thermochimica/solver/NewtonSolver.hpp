/// @file NewtonSolver.hpp
/// @brief Newton direction solver with reduced active element system
/// @details Implements INewtonSolver using reduced Hessian with active elements only

#pragma once

#include "thermochimica/interfaces/INewtonSolver.hpp"
#include <Eigen/Dense>
#include <vector>

namespace Thermochimica {

// Forward declarations
class ThermoState;
class GEMState;
class ThermoIO;
class PhaseConstraints;

/// @brief Newton solver with reduced system (active elements only)
/// @details Constructs Hessian and RHS for reduced Newton system excluding
/// elements with zero input mass. Optionally includes phase constraint contributions.
class NewtonSolver : public INewtonSolver {
public:
    /// @brief Constructor
    /// @param io Reference to ThermoIO for temperature
    /// @param constraints Pointer to PhaseConstraints (can be nullptr)
    NewtonSolver(ThermoIO& io, PhaseConstraints* constraints = nullptr);

    /// @brief Compute Newton direction using reduced system
    /// @param state Thermodynamic state
    /// @param gemState GEM solver state
    /// @param direction Output: Newton direction vector
    /// @return Error code (0 = success, 1 = NaN/Inf detected)
    int computeDirection(ThermoState& state,
                        GEMState& gemState,
                        Eigen::VectorXd& direction) override;

    /// @brief Handle singular matrix by adding diagonal regularization
    /// @param hessian Hessian matrix (modified in-place)
    /// @param rhs RHS vector (modified in-place)
    void handleSingular(Eigen::MatrixXd& hessian,
                       Eigen::VectorXd& rhs) override;

    /// @brief Get solver name
    /// @return "NewtonSolver"
    const char* getNewtonSolverName() const override {
        return "NewtonSolver";
    }

private:
    /// @brief Build list of active elements (non-zero input mass)
    /// @param state Thermodynamic state
    /// @param activeElements Output: indices of active elements
    void buildActiveElements(const ThermoState& state,
                            std::vector<int>& activeElements) const;

    /// @brief Construct reduced Hessian and RHS
    /// @param state Thermodynamic state
    /// @param hessian Output: reduced Hessian matrix
    /// @param rhs Output: reduced RHS vector
    /// @param activeElements List of active element indices
    void constructReducedSystem(const ThermoState& state,
                               Eigen::MatrixXd& hessian,
                               Eigen::VectorXd& rhs,
                               const std::vector<int>& activeElements) const;

    /// @brief Add phase constraint contributions to Newton system
    /// @param state Thermodynamic state
    /// @param gemState GEM solver state
    /// @param hessian Hessian matrix (modified in-place)
    /// @param rhs RHS vector (modified in-place)
    /// @param activeElements List of active element indices
    void addConstraintContributions(const ThermoState& state,
                                   const GEMState& gemState,
                                   Eigen::MatrixXd& hessian,
                                   Eigen::VectorXd& rhs,
                                   const std::vector<int>& activeElements) const;

    /// @brief Map reduced solution back to full update vector
    /// @param updateReduced Reduced solution vector
    /// @param direction Output: full direction vector
    /// @param state Thermodynamic state
    /// @param activeElements List of active element indices
    void mapToFullDirection(const Eigen::VectorXd& updateReduced,
                           Eigen::VectorXd& direction,
                           const ThermoState& state,
                           const std::vector<int>& activeElements) const;

    ThermoIO& io_;                    ///< Reference to ThermoIO
    PhaseConstraints* constraints_;   ///< Pointer to PhaseConstraints (optional)
};

} // namespace Thermochimica
