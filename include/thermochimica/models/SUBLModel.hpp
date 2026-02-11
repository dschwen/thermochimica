/// @file SUBLModel.hpp
/// @brief Compound Energy Formalism (SUBL/SUBLM) model implementing IThermodynamicModel
/// @details CEF model for phases with sublattice structure

#pragma once

#include "thermochimica/interfaces/IThermodynamicModel.hpp"

namespace Thermochimica {

// Forward declarations are in IThermodynamicModel.hpp

/// @brief SUBL (Compound Energy Formalism) thermodynamic model
/// @details Implements IThermodynamicModel for SUBL and SUBLM phase types.
/// Uses bridge pattern to call legacy implementation.
///
/// The Compound Energy Formalism (CEF) treats phases with sublattice structure,
/// commonly used for intermetallic compounds and ionic phases.
class SUBLModel : public IThermodynamicModel {
public:
    /// @brief Constructor
    SUBLModel() = default;

    /// @brief Destructor
    ~SUBLModel() override = default;

    /// @brief Compute excess Gibbs energy for SUBL phase
    /// @param state Thermodynamic state
    /// @param gemState GEM solver state
    /// @param phaseIndex Phase index (0-based)
    void computeExcessGibbs(ThermoState& state, GEMState& gemState, int phaseIndex) override;

    /// @brief Get model type
    /// @return PhaseType::SUBL
    Constants::PhaseType getModelType() const override {
        return Constants::PhaseType::SUBL;
    }

    /// @brief Get model name
    /// @return "SUBLModel"
    const char* getModelName() const override {
        return "SUBLModel";
    }

    /// @brief Check if this model can handle the given phase
    /// @param state Thermodynamic state
    /// @param phaseIndex Phase index (0-based)
    /// @return true if phase is SUBL or SUBLM type
    bool canHandle(const ThermoState& state, int phaseIndex) const override;
};

} // namespace Thermochimica
