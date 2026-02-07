/// @file WolfeLineSearch.hpp
/// @brief Wolfe conditions line search strategy
/// @details Implements ILineSearch using Wolfe conditions for step length selection

#pragma once

#include "thermochimica/interfaces/ILineSearch.hpp"
#include <Eigen/Dense>

namespace Thermochimica {

// Forward declarations
class ThermoState;
class GEMState;

/// @brief Wolfe conditions line search
/// @details Backtracking line search with Wolfe conditions (sufficient decrease + curvature)
class WolfeLineSearch : public ILineSearch {
public:
    /// @brief Perform line search along Newton direction
    /// @param state Thermodynamic state
    /// @param gemState GEM solver state
    /// @param direction Newton direction vector
    /// @param stepLength Output: computed step length
    void search(ThermoState& state,
               GEMState& gemState,
               const Eigen::VectorXd& direction,
               double& stepLength) override;

    /// @brief Update state variables with computed step
    /// @param state Thermodynamic state to update
    /// @param gemState GEM solver state to update
    /// @param direction Newton direction vector
    /// @param stepLength Step length
    void updateState(ThermoState& state,
                    GEMState& gemState,
                    const Eigen::VectorXd& direction,
                    double stepLength) override;

    /// @brief Get line search name
    /// @return "WolfeLineSearch"
    const char* getLineSearchName() const override {
        return "WolfeLineSearch";
    }

private:
    /// @brief Initialize step length
    void initStepLength(ThermoState& state, GEMState& gemState, double& stepLength);

    /// @brief Update element potentials
    void updateElementPotentials(ThermoState& state, GEMState& gemState,
                                 const Eigen::VectorXd& direction, double stepLength);

    /// @brief Update phase moles
    void updatePhaseMoles(ThermoState& state, GEMState& gemState,
                         const Eigen::VectorXd& direction, double stepLength);

    /// @brief Update species moles and fractions
    void updateSpecies(ThermoState& state);

    /// @brief Check Wolfe conditions
    bool checkWolfeConditions(const GEMState& gemState, double functionNormLast) const;
};

} // namespace Thermochimica
