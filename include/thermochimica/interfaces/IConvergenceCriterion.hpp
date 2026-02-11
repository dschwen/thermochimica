/// @file IConvergenceCriterion.hpp
/// @brief Interface for individual convergence criteria
/// @details Defines the contract for convergence checking components

#pragma once

#include <string>

namespace Thermochimica {

// Forward declarations
class ThermoState;
class GEMState;

/// @brief Abstract interface for convergence criteria
/// @details Implementations check specific convergence conditions:
/// - MassBalanceCheck: Element mass balance satisfaction
/// - ChemicalPotentialCheck: Chemical potential equilibrium
/// - PhaseRuleCheck: Phase rule compliance
/// - DrivingForceCheck: Unstable phase driving forces
/// - MiscibilityGapCheck: Miscibility gap detection
/// - ConstraintCheck: Phase fraction constraint satisfaction
class IConvergenceCriterion {
public:
    virtual ~IConvergenceCriterion() = default;

    /// @brief Check if this convergence criterion is satisfied
    /// @param state Thermodynamic state
    /// @param gemState GEM solver state
    /// @return true if criterion is satisfied
    virtual bool check(const ThermoState& state,
                      const GEMState& gemState) const = 0;

    /// @brief Get criterion name for logging
    /// @return Criterion name string (e.g., "MassBalance", "ChemicalPotential")
    virtual std::string getName() const = 0;

    /// @brief Get detailed status message
    /// @param state Thermodynamic state
    /// @param gemState GEM solver state
    /// @return Status message describing current state of criterion
    virtual std::string getStatusMessage(const ThermoState& state,
                                        const GEMState& gemState) const = 0;

    /// @brief Get tolerance value for this criterion
    /// @return Tolerance value (criterion-specific)
    virtual double getTolerance() const = 0;
};

} // namespace Thermochimica
