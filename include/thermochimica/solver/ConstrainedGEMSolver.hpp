/// @file ConstrainedGEMSolver.hpp
/// @brief Constrained GEM solver with phase fraction targets
/// @details Implements ISolver for phase-constrained Gibbs energy minimization
/// using augmented Lagrangian method for phase field modeling

#pragma once

#include "thermochimica/interfaces/ISolver.hpp"

namespace Thermochimica {

// Forward declarations
class ThermoState;
class ThermoIO;
class GEMState;
class PhaseConstraints;

/// @brief Constrained GEM solver with phase fraction targets
/// @details Uses augmented Lagrangian method to enforce phase fraction constraints
/// while minimizing Gibbs energy. Suitable for phase field modeling where specific
/// phase fractions must be maintained. Currently bridges to legacy implementation.
class ConstrainedGEMSolver : public ISolver {
public:
    /// @brief Constructor
    /// @param constraints Reference to phase constraints (must outlive solver)
    explicit ConstrainedGEMSolver(PhaseConstraints& constraints);

    /// @brief Run constrained solver to convergence
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
    /// @details Checks both inner loop convergence AND constraint satisfaction
    /// @param state Thermodynamic state
    /// @param gemState GEM solver state
    /// @return true if converged and constraints satisfied
    bool isConverged(const ThermoState& state,
                    const GEMState& gemState) const override;

    /// @brief Get solver name
    /// @return "ConstrainedGEMSolver"
    const char* getSolverName() const override {
        return "ConstrainedGEMSolver";
    }

private:
    PhaseConstraints& constraints_;  ///< Reference to phase constraints
};

} // namespace Thermochimica
