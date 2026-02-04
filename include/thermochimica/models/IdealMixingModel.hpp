/// @file IdealMixingModel.hpp
/// @brief Ideal mixing thermodynamic model
/// @details Implements IThermodynamicModel for ideal solution phases (IDMX)

#pragma once

#include "thermochimica/interfaces/IThermodynamicModel.hpp"

namespace Thermochimica {

// Forward declarations
class ThermoState;
class GEMState;

/// @brief Ideal mixing model (no excess Gibbs energy)
/// @details For ideal solutions: G^xs = 0, only configurational entropy RT*ln(x)
class IdealMixingModel : public IThermodynamicModel {
public:
    /// @brief Compute excess Gibbs energy for ideal mixing
    /// @param state Thermodynamic state
    /// @param gemState GEM solver state
    /// @param phaseIndex Phase index
    void computeExcessGibbs(ThermoState& state,
                           GEMState& gemState,
                           int phaseIndex) override;

    /// @brief Get model type
    /// @return PhaseType::IDMX
    Constants::PhaseType getModelType() const override {
        return Constants::PhaseType::IDMX;
    }

    /// @brief Check if model can handle this phase
    /// @param state Thermodynamic state
    /// @param phaseIndex Phase index
    /// @return true if phase type is IDMX
    bool canHandle(const ThermoState& state, int phaseIndex) const override;

    /// @brief Get model name
    /// @return "IdealMixing"
    const char* getModelName() const override {
        return "IdealMixing";
    }
};

} // namespace Thermochimica
