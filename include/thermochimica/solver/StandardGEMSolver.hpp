/// @file StandardGEMSolver.hpp
/// @brief Standard unconstrained GEM solver implementation
/// @details Implements ISolver for unconstrained Gibbs energy minimization

#pragma once

#include "thermochimica/interfaces/ISolver.hpp"
#include <memory>

namespace Thermochimica {

// Forward declarations
class ThermoState;
class ThermoIO;
class GEMState;
class PhaseConstraints;

/// @brief Standard unconstrained GEM solver
/// @details Minimizes Gibbs energy without phase fraction constraints
/// using Newton's method with line search and phase assemblage management.
/// Currently bridges to legacy GEMSolver implementation.
class StandardGEMSolver : public ISolver {
public:
    /// @brief Constructor
    StandardGEMSolver();

    /// @brief Run solver to convergence
    /// @param state Thermodynamic state
    /// @param io Input/output parameters
    /// @param gemState GEM solver iteration state
    /// @param phaseManager Phase assemblage manager
    /// @param newton Newton solver for direction computation
    /// @param lineSearch Line search strategy
    /// @param models Vector of thermodynamic models
    /// @return Error code (0 = success)
    int solve(ThermoState& state,
             ThermoIO& io,
             GEMState& gemState,
             PhaseAssemblageManager& phaseManager,
             INewtonSolver& newton,
             ILineSearch& lineSearch,
             const std::vector<IThermodynamicModel*>& models) override;

    /// @brief Initialize solver state before iteration
    /// @param state Thermodynamic state
    /// @param gemState GEM solver state
    void initialize(ThermoState& state, GEMState& gemState) override;

    /// @brief Check if solver has converged
    /// @param state Thermodynamic state
    /// @param gemState GEM solver state
    /// @return true if converged
    bool isConverged(const ThermoState& state,
                    const GEMState& gemState) const override;

    /// @brief Get solver name
    /// @return "StandardGEMSolver"
    const char* getSolverName() const override {
        return "StandardGEMSolver";
    }

private:
    std::unique_ptr<PhaseConstraints> phaseConstraints_;  ///< Owned PhaseConstraints for unconstrained solve
};

} // namespace Thermochimica
