/// @file SUBQModel.hpp
/// @brief Modified Quasichemical Model in Quadruplet (SUBQ) implementing IThermodynamicModel
/// @details MQM quadruplet approximation for molten systems

#pragma once

#include "thermochimica/interfaces/IThermodynamicModel.hpp"

namespace Thermochimica {

// Forward declarations are in IThermodynamicModel.hpp

/// @brief SUBQ (Modified Quasichemical Model - Quadruplet) thermodynamic model
/// @details Implements IThermodynamicModel for SUBQ phase type.
/// Uses bridge pattern to call legacy implementation.
///
/// The Modified Quasichemical Model in the Quadruplet approximation
/// is similar to SUBG but uses different coordination assumptions.
class SUBQModel : public IThermodynamicModel {
public:
    /// @brief Constructor
    SUBQModel() = default;

    /// @brief Destructor
    ~SUBQModel() override = default;

    /// @brief Compute excess Gibbs energy for SUBQ phase
    /// @param state Thermodynamic state
    /// @param gemState GEM solver state
    /// @param phaseIndex Phase index (0-based)
    void computeExcessGibbs(ThermoState& state, GEMState& gemState, int phaseIndex) override;

    /// @brief Get model type
    /// @return PhaseType::SUBQ
    Constants::PhaseType getModelType() const override {
        return Constants::PhaseType::SUBQ;
    }

    /// @brief Get model name
    /// @return "SUBQModel"
    const char* getModelName() const override {
        return "SUBQModel";
    }

    /// @brief Check if this model can handle the given phase
    /// @param state Thermodynamic state
    /// @param phaseIndex Phase index (0-based)
    /// @return true if phase is SUBQ type
    bool canHandle(const ThermoState& state, int phaseIndex) const override;
};

} // namespace Thermochimica
