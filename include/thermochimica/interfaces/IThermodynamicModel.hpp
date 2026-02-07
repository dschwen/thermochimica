/// @file IThermodynamicModel.hpp
/// @brief Interface for thermodynamic models (QKTO, RKMP, SUBL, etc.)
/// @details Defines the contract for computing excess Gibbs energy contributions

#pragma once

#include "thermochimica/util/Constants.hpp"

namespace Thermochimica {

// Forward declarations
class ThermoState;
class GEMState;

/// @brief Abstract interface for thermodynamic models
/// @details Implementations compute excess Gibbs energy for solution phases
/// using different mixing models (ideal, QKTO, RKMP, SUBL, SUBG, SUBQ, etc.)
class IThermodynamicModel {
public:
    virtual ~IThermodynamicModel() = default;

    /// @brief Compute excess Gibbs energy for a phase
    /// @param state Thermodynamic state containing species data
    /// @param gemState GEM solver state for storing partial excess Gibbs
    /// @param phaseIndex Index of the phase to compute
    virtual void computeExcessGibbs(ThermoState& state,
                                    GEMState& gemState,
                                    int phaseIndex) = 0;

    /// @brief Get the model type identifier
    /// @return PhaseType enum value
    virtual Constants::PhaseType getModelType() const = 0;

    /// @brief Check if this model can handle the given phase
    /// @param state Thermodynamic state
    /// @param phaseIndex Phase to check
    /// @return true if this model applies to the phase
    virtual bool canHandle(const ThermoState& state, int phaseIndex) const = 0;

    /// @brief Get human-readable model name
    /// @return Model name string
    virtual const char* getModelName() const = 0;
};

} // namespace Thermochimica
