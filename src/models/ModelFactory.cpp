/// @file ModelFactory.cpp
/// @brief Implementation of thermodynamic model factory

#include "thermochimica/models/ModelFactory.hpp"
#include "thermochimica/models/IdealMixingModel.hpp"
#include "thermochimica/models/QKTOModel.hpp"
// TODO: Add other model headers when implemented
// #include "thermochimica/models/RKMPModel.hpp"
// #include "thermochimica/models/SUBLModel.hpp"
// #include "thermochimica/models/SUBGModel.hpp"
// #include "thermochimica/models/SUBQModel.hpp"

namespace Thermochimica {

ModelFactory::ModelFactory() {
    registerStandardModels();
}

IThermodynamicModel* ModelFactory::getModel(Constants::PhaseType phaseType) {
    auto it = models_.find(phaseType);
    if (it != models_.end()) {
        return it->second.get();
    }
    return nullptr;
}

void ModelFactory::registerModel(Constants::PhaseType phaseType,
                                 std::unique_ptr<IThermodynamicModel> model) {
    models_[phaseType] = std::move(model);
}

bool ModelFactory::hasModel(Constants::PhaseType phaseType) const {
    return models_.find(phaseType) != models_.end();
}

std::vector<Constants::PhaseType> ModelFactory::getRegisteredTypes() const {
    std::vector<Constants::PhaseType> types;
    types.reserve(models_.size());
    for (const auto& pair : models_) {
        types.push_back(pair.first);
    }
    return types;
}

void ModelFactory::registerStandardModels() {
    // Register ideal mixing model
    registerModel(Constants::PhaseType::IDMX,
                 std::make_unique<IdealMixingModel>());

    // Register QKTO model
    registerModel(Constants::PhaseType::QKTO,
                 std::make_unique<QKTOModel>());

    // TODO: Register other models when implemented
    // registerModel(Constants::PhaseType::RKMP, std::make_unique<RKMPModel>());
    // registerModel(Constants::PhaseType::RKMPM, std::make_unique<RKMPModel>());
    // registerModel(Constants::PhaseType::SUBL, std::make_unique<SUBLModel>());
    // registerModel(Constants::PhaseType::SUBLM, std::make_unique<SUBLModel>());
    // registerModel(Constants::PhaseType::SUBG, std::make_unique<SUBGModel>());
    // registerModel(Constants::PhaseType::SUBQ, std::make_unique<SUBQModel>());
}

} // namespace Thermochimica
