/// @file SUBGModel.hpp
/// @brief Modified Quasichemical Model (SUBG) implementing IThermodynamicModel
/// @details MQM model for phases with quadruplet interactions

#pragma once

#include "thermochimica/interfaces/IThermodynamicModel.hpp"

namespace Thermochimica {

// Forward declarations are in IThermodynamicModel.hpp

/// @brief SUBG (Modified Quasichemical Model) thermodynamic model
/// @details Implements IThermodynamicModel for SUBG phase type.
/// Uses bridge pattern to call legacy implementation.
///
/// The Modified Quasichemical Model (MQM) in the quadruplet approximation
/// treats nearest-neighbor interactions for molten salts, oxides, and slags.
class SUBGModel : public IThermodynamicModel {
public:
    /// @brief Constructor
    SUBGModel() = default;

    /// @brief Destructor
    ~SUBGModel() override = default;

    /// @brief Compute excess Gibbs energy for SUBG phase
    /// @param state Thermodynamic state
    /// @param gemState GEM solver state
    /// @param phaseIndex Phase index (0-based)
    void computeExcessGibbs(ThermoState& state, GEMState& gemState, int phaseIndex) override;

    /// @brief Get model type
    /// @return PhaseType::SUBG
    Constants::PhaseType getModelType() const override {
        return Constants::PhaseType::SUBG;
    }

    /// @brief Get model name
    /// @return "SUBGModel"
    const char* getModelName() const override {
        return "SUBGModel";
    }

    /// @brief Check if this model can handle the given phase
    /// @param state Thermodynamic state
    /// @param phaseIndex Phase index (0-based)
    /// @return true if phase is SUBG type
    bool canHandle(const ThermoState& state, int phaseIndex) const override;
};

} // namespace Thermochimica
