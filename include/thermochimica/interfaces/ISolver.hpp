/// @file ISolver.hpp
/// @brief Interface for GEM solver strategies
/// @details Defines the contract for different solver implementations
/// (standard GEM, constrained GEM, etc.)

#pragma once

#include <vector>

namespace Thermochimica {

// Forward declarations
class ThermoState;
class ThermoIO;
class GEMState;
class PhaseAssemblageManager;
class INewtonSolver;
class ILineSearch;
class IThermodynamicModel;

/// @brief Abstract interface for GEM solver strategies
/// @details Implementations provide different solving approaches:
/// - StandardGEMSolver: Unconstrained Gibbs energy minimization
/// - ConstrainedGEMSolver: Constrained minimization with phase fraction targets
class ISolver {
public:
    virtual ~ISolver() = default;

    /// @brief Run solver to convergence
    /// @param state Thermodynamic state
    /// @param io Input/output parameters
    /// @param gemState GEM solver iteration state
    /// @param phaseManager Phase assemblage manager
    /// @param newton Newton solver for direction computation
    /// @param lineSearch Line search strategy
    /// @param models Vector of thermodynamic models
    /// @return Error code (0 = success)
    virtual int solve(ThermoState& state,
                     ThermoIO& io,
                     GEMState& gemState,
                     PhaseAssemblageManager& phaseManager,
                     INewtonSolver& newton,
                     ILineSearch& lineSearch,
                     const std::vector<IThermodynamicModel*>& models) = 0;

    /// @brief Initialize solver state before iteration
    /// @param state Thermodynamic state
    /// @param gemState GEM solver state
    virtual void initialize(ThermoState& state, GEMState& gemState) = 0;

    /// @brief Check if solver has converged
    /// @param state Thermodynamic state
    /// @param gemState GEM solver state
    /// @return true if converged
    virtual bool isConverged(const ThermoState& state,
                            const GEMState& gemState) const = 0;

    /// @brief Get solver name for logging
    /// @return Solver name string
    virtual const char* getSolverName() const = 0;
};

} // namespace Thermochimica
