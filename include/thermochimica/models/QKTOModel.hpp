/// @file QKTOModel.hpp
/// @brief Quasi-chemical Kohler-Toop (QKTO) thermodynamic model
/// @details Implements IThermodynamicModel for QKTO solution phases

#pragma once

#include "thermochimica/interfaces/IThermodynamicModel.hpp"

namespace Thermochimica {

// Forward declarations
class ThermoState;
class GEMState;
class ThermoContext;

/// @brief Quasi-chemical Kohler-Toop model
/// @details Polynomial regular solution with Kohler interpolation for ternary systems
class QKTOModel : public IThermodynamicModel {
public:
    /// @brief Compute excess Gibbs energy for QKTO phase
    /// @param state Thermodynamic state
    /// @param gemState GEM solver state
    /// @param phaseIndex Phase index
    void computeExcessGibbs(ThermoState& state,
                           GEMState& gemState,
                           int phaseIndex) override;

    /// @brief Get model type
    /// @return PhaseType::QKTO
    Constants::PhaseType getModelType() const override {
        return Constants::PhaseType::QKTO;
    }

    /// @brief Check if model can handle this phase
    /// @param state Thermodynamic state
    /// @param phaseIndex Phase index
    /// @return true if phase type is QKTO
    bool canHandle(const ThermoState& state, int phaseIndex) const override;

    /// @brief Get model name
    /// @return "QKTO"
    const char* getModelName() const override {
        return "QKTO";
    }
};

} // namespace Thermochimica
