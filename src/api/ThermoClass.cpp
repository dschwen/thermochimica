/// @file ThermoClass.cpp
/// @brief Implementation of main Thermochimica class

#include "thermochimica/ThermoClass.hpp"
#include "thermochimica/Thermochimica.hpp"  // For legacy free functions

// Interface headers (needed for unique_ptr destruction)
#include "thermochimica/interfaces/ISolver.hpp"
#include "thermochimica/interfaces/INewtonSolver.hpp"
#include "thermochimica/interfaces/ILineSearch.hpp"
#include "thermochimica/interfaces/IParser.hpp"
#include "thermochimica/interfaces/IThermodynamicModel.hpp"

// Implementation headers
#include "thermochimica/solver/StandardGEMSolver.hpp"
#include "thermochimica/solver/ConstrainedGEMSolver.hpp"
#include "thermochimica/solver/NewtonSolver.hpp"
#include "thermochimica/solver/WolfeLineSearch.hpp"
#include "thermochimica/solver/PhaseAssemblageManager.hpp"
#include "thermochimica/models/ModelFactory.hpp"
#include "thermochimica/context/ThermoState.hpp"
#include "thermochimica/context/ThermoIO.hpp"
#include "thermochimica/context/GEMState.hpp"
#include "thermochimica/context/PhaseConstraints.hpp"

#include <vector>

namespace Thermochimica {

ThermoClass::ThermoClass()
    : state_(context_.thermo.get()),
      io_(context_.io.get()),
      gemState_(context_.gem.get()),
      phaseConstraints_(context_.phaseConstraints.get()) {

    initializeDefaultStrategies();
}

ThermoClass::~ThermoClass() = default;

ThermoClass::ThermoClass(ThermoClass&&) noexcept = default;

ThermoClass& ThermoClass::operator=(ThermoClass&&) noexcept = default;

void ThermoClass::initializeDefaultStrategies() {
    // Initialize model factory
    modelFactory_ = std::make_unique<ModelFactory>();

    // Initialize default Newton solver
    newtonSolver_ = std::make_unique<NewtonSolver>(*io_, phaseConstraints_);

    // Initialize default line search
    lineSearch_ = std::make_unique<WolfeLineSearch>();

    // Initialize phase assemblage manager
    phaseManager_ = std::make_unique<PhaseAssemblageManager>(*state_, *gemState_, *io_);

    // Initialize default solver (standard unconstrained)
    solver_ = std::make_unique<StandardGEMSolver>();
}

// =========================================================================
// Database Loading
// =========================================================================

void ThermoClass::setThermoFilename(const std::string& filename) {
    Thermochimica::setThermoFilename(context_, filename);
}

int ThermoClass::parseCSDataFile() {
    Thermochimica::parseCSDataFile(context_);
    return context_.infoThermo();
}

int ThermoClass::loadDatabase(const std::string& filename) {
    setThermoFilename(filename);
    parseCSDataFile();
    return context_.infoThermo();
}

void ThermoClass::setParser(std::unique_ptr<IParser> parser) {
    parser_ = std::move(parser);
}

// =========================================================================
// Input Configuration
// =========================================================================

void ThermoClass::setTemperature(double temperature) {
    Thermochimica::setTemperature(context_, temperature);
}

void ThermoClass::setPressure(double pressure) {
    Thermochimica::setPressure(context_, pressure);
}

void ThermoClass::setTemperaturePressure(double temperature, double pressure) {
    Thermochimica::setTemperaturePressure(context_, temperature, pressure);
}

void ThermoClass::setElementMass(int atomicNumber, double mass) {
    Thermochimica::setElementMass(context_, atomicNumber, mass);
}

void ThermoClass::setElementMass(const std::string& elementName, double mass) {
    Thermochimica::setElementMass(context_, elementName, mass);
}

void ThermoClass::setStandardUnits() {
    Thermochimica::setStandardUnits(context_);
}

void ThermoClass::setUnitsSI() {
    // Set SI units: K, Pa, moles
    Thermochimica::setUnits(context_, "K", "Pa", "moles");
}

// =========================================================================
// Solver Configuration
// =========================================================================

void ThermoClass::setSolver(std::unique_ptr<ISolver> solver) {
    solver_ = std::move(solver);
}

void ThermoClass::setNewtonSolver(std::unique_ptr<INewtonSolver> newton) {
    newtonSolver_ = std::move(newton);
}

void ThermoClass::setLineSearch(std::unique_ptr<ILineSearch> lineSearch) {
    lineSearch_ = std::move(lineSearch);
}

// =========================================================================
// Phase Constraints
// =========================================================================

void ThermoClass::setSolnPhaseConstraint(const std::string& phaseName, double targetFraction) {
    Thermochimica::setSolnPhaseConstraint(context_, phaseName, targetFraction);

    // Switch to constrained solver if constraints are active
    if (phaseConstraints_->hasActiveConstraints()) {
        solver_ = std::make_unique<ConstrainedGEMSolver>(*phaseConstraints_);
    }
}

void ThermoClass::setCondPhaseConstraint(const std::string& speciesName, double targetFraction) {
    Thermochimica::setCondPhaseConstraint(context_, speciesName, targetFraction);

    // Switch to constrained solver if constraints are active
    if (phaseConstraints_->hasActiveConstraints()) {
        solver_ = std::make_unique<ConstrainedGEMSolver>(*phaseConstraints_);
    }
}

void ThermoClass::removePhaseConstraint(const std::string& phaseName) {
    Thermochimica::removePhaseConstraint(context_, phaseName);

    // Switch back to standard solver if no constraints remain
    if (!phaseConstraints_->hasActiveConstraints()) {
        solver_ = std::make_unique<StandardGEMSolver>();
    }
}

void ThermoClass::clearPhaseConstraints() {
    Thermochimica::clearPhaseConstraints(context_);

    // Switch back to standard solver
    solver_ = std::make_unique<StandardGEMSolver>();
}

// =========================================================================
// Main Computation
// =========================================================================

int ThermoClass::calculate() {
    initialize();
    if (context_.infoThermo() != 0) {
        return context_.infoThermo();
    }

    checkSystem();
    if (context_.infoThermo() != 0) {
        return context_.infoThermo();
    }

    computeThermoData();
    if (context_.infoThermo() != 0) {
        return context_.infoThermo();
    }

    setup();
    if (context_.infoThermo() != 0) {
        return context_.infoThermo();
    }

    return solve();
}

// =========================================================================
// Granular Control
// =========================================================================

void ThermoClass::initialize() {
    Thermochimica::init(context_);
}

void ThermoClass::checkSystem() {
    Thermochimica::checkSystem(context_);
}

void ThermoClass::computeThermoData() {
    Thermochimica::compThermoData(context_);
}

void ThermoClass::setup() {
    Thermochimica::setup(context_);
}

int ThermoClass::solve() {
    if (!solver_) {
        solver_ = std::make_unique<StandardGEMSolver>();
    }
    if (!newtonSolver_) {
        newtonSolver_ = std::make_unique<NewtonSolver>(*io_, phaseConstraints_);
    }
    if (!lineSearch_) {
        lineSearch_ = std::make_unique<WolfeLineSearch>();
    }
    if (!phaseManager_) {
        phaseManager_ = std::make_unique<PhaseAssemblageManager>(*state_, *gemState_, *io_);
    }
    if (!modelFactory_) {
        modelFactory_ = std::make_unique<ModelFactory>();
    }

    std::vector<IThermodynamicModel*> models;
    for (const auto type : modelFactory_->getRegisteredTypes()) {
        if (auto* model = modelFactory_->getModel(type)) {
            models.push_back(model);
        }
    }

    int result = solver_->solve(*state_, *io_, *gemState_, *phaseManager_,
                                *newtonSolver_, *lineSearch_, models);

    if (result != 0) {
        context_.setInfoThermo(result);
        return context_.infoThermo();
    }

    Thermochimica::postProcess(context_);
    return context_.infoThermo();
}

// =========================================================================
// Output Retrieval
// =========================================================================

std::pair<double, int> ThermoClass::getOutputChemPot(const std::string& elementName) const {
    return Thermochimica::getOutputChemPot(context_, elementName);
}

std::pair<double, int> ThermoClass::getMolesPhase(const std::string& phaseName) const {
    return Thermochimica::getMolesPhase(context_, phaseName);
}

double ThermoClass::getGibbsEnergy() const {
    return Thermochimica::getGibbsEnergy(context_);
}

void ThermoClass::printResults() const {
    // printResults modifies context, so we need const_cast
    Thermochimica::printResults(const_cast<ThermoContext&>(context_));
}

// =========================================================================
// Status
// =========================================================================

int ThermoClass::getInfoCode() const {
    return context_.infoThermo();
}

bool ThermoClass::isSuccess() const {
    return context_.infoThermo() == 0;
}

std::string ThermoClass::getErrorMessage() const {
    return Thermochimica::getErrorMessage(context_.infoThermo());
}

// =========================================================================
// Reset
// =========================================================================

void ThermoClass::reset() {
    Thermochimica::resetThermo(context_);
}

void ThermoClass::resetAll() {
    Thermochimica::resetThermoAll(context_);
    initializeDefaultStrategies();
}

void ThermoClass::validateConfiguration() const {
    // Basic validation - could be expanded
    if (io_->dTemperature <= 0) {
        throw std::runtime_error("Temperature not set or invalid");
    }
    if (io_->dPressure <= 0) {
        throw std::runtime_error("Pressure not set or invalid");
    }
}

} // namespace Thermochimica
