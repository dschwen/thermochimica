/// @file ModelFactory.hpp
/// @brief Factory for creating thermodynamic models
/// @details Manages creation and lifetime of IThermodynamicModel implementations

#pragma once

#include "thermochimica/interfaces/IThermodynamicModel.hpp"
#include "thermochimica/util/Constants.hpp"
#include <map>
#include <memory>

namespace Thermochimica {

/// @brief Factory for thermodynamic models
/// @details Creates and manages model instances based on PhaseType.
/// Supports registration of custom models at runtime.
class ModelFactory {
public:
    /// @brief Constructor - registers all standard models
    ModelFactory();

    /// @brief Destructor
    ~ModelFactory() = default;

    /// @brief Get model for a phase type
    /// @param phaseType Type of phase
    /// @return Pointer to model (nullptr if not found)
    IThermodynamicModel* getModel(Constants::PhaseType phaseType);

    /// @brief Register a custom model
    /// @param phaseType Type of phase this model handles
    /// @param model Unique pointer to model (factory takes ownership)
    void registerModel(Constants::PhaseType phaseType,
                      std::unique_ptr<IThermodynamicModel> model);

    /// @brief Check if factory has a model for this phase type
    /// @param phaseType Type of phase
    /// @return true if model is registered
    bool hasModel(Constants::PhaseType phaseType) const;

    /// @brief Get all registered phase types
    /// @return Vector of phase types with registered models
    std::vector<Constants::PhaseType> getRegisteredTypes() const;

private:
    /// Map from phase type to model instance
    std::map<Constants::PhaseType, std::unique_ptr<IThermodynamicModel>> models_;

    /// Register all standard models (called by constructor)
    void registerStandardModels();
};

} // namespace Thermochimica
