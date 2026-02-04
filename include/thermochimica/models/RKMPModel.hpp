/// @file RKMPModel.hpp
/// @brief Redlich-Kister-Muggianu-Polynomial model implementing IThermodynamicModel
/// @details RKMP/RKMPM model for non-ideal solutions with binary, ternary, and quaternary interactions

#pragma once

#include "thermochimica/interfaces/IThermodynamicModel.hpp"

namespace Thermochimica {

// Forward declarations are in IThermodynamicModel.hpp

/// @brief RKMP (Redlich-Kister-Muggianu-Polynomial) thermodynamic model
/// @details Implements IThermodynamicModel for RKMP and RKMPM phase types.
/// Uses bridge pattern to call legacy implementation.
///
/// Binary interactions:  g^ex = x1*x2 * sum L_{1,2} * (x1-x2)^v
/// Ternary interactions: g^ex = x1*x2*x3 * L_{1,2,3}
/// Quaternary: g^ex = x1*x2*x3*x4 * L_{1,2,3,4}
class RKMPModel : public IThermodynamicModel {
public:
    /// @brief Constructor
    RKMPModel() = default;

    /// @brief Destructor
    ~RKMPModel() override = default;

    /// @brief Compute excess Gibbs energy for RKMP phase
    /// @param state Thermodynamic state
    /// @param gemState GEM solver state
    /// @param phaseIndex Phase index (0-based)
    void computeExcessGibbs(ThermoState& state, GEMState& gemState, int phaseIndex) override;

    /// @brief Get model type
    /// @return PhaseType::RKMP
    Constants::PhaseType getModelType() const override {
        return Constants::PhaseType::RKMP;
    }

    /// @brief Get model name
    /// @return "RKMPModel"
    const char* getModelName() const override {
        return "RKMPModel";
    }

    /// @brief Check if this model can handle the given phase
    /// @param state Thermodynamic state
    /// @param phaseIndex Phase index (0-based)
    /// @return true if phase is RKMP or RKMPM type
    bool canHandle(const ThermoState& state, int phaseIndex) const override;
};

} // namespace Thermochimica
