/// @file ThermoClass.cpp
/// @brief Implementation of main Thermochimica class

#include "thermochimica/ThermoClass.hpp"

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
#include "thermochimica/parser/ChemSageParser.hpp"
#include "thermochimica/util/ErrorCodes.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace Thermochimica {

// Forward declarations for internal functions (declared in other compilation units)
void levelingSolver(ThermoContext& ctx);
void postLevelingSolver(ThermoContext& ctx);
void postProcess(ThermoContext& ctx);

ThermoClass::ThermoClass()
    : state_(context_.thermo.get()),
      io_(context_.io.get()),
      gemState_(context_.gem.get()),
      phaseConstraints_(context_.phaseConstraints.get()) {

    initializeDefaultStrategies();
}

ThermoClass::~ThermoClass() = default;

ThermoClass::ThermoClass(ThermoClass&& other) noexcept
    : context_(std::move(other.context_)),
      state_(context_.thermo.get()),
      io_(context_.io.get()),
      gemState_(context_.gem.get()),
      phaseConstraints_(context_.phaseConstraints.get()),
      parser_(std::move(other.parser_)),
      solver_(std::move(other.solver_)),
      newtonSolver_(std::move(other.newtonSolver_)),
      lineSearch_(std::move(other.lineSearch_)),
      phaseManager_(std::move(other.phaseManager_)),
      modelFactory_(std::move(other.modelFactory_)) {

    // Rebuild strategies that hold references to state objects
    // (they were constructed with references to the old context)
    if (newtonSolver_) {
        newtonSolver_ = std::make_unique<NewtonSolver>(*io_, phaseConstraints_);
    }
    if (phaseManager_) {
        phaseManager_ = std::make_unique<PhaseAssemblageManager>(*state_, *gemState_, *io_);
    }
}

ThermoClass& ThermoClass::operator=(ThermoClass&& other) noexcept {
    if (this != &other) {
        // Move the context and unique_ptr members
        context_ = std::move(other.context_);
        parser_ = std::move(other.parser_);
        solver_ = std::move(other.solver_);
        newtonSolver_ = std::move(other.newtonSolver_);
        lineSearch_ = std::move(other.lineSearch_);
        phaseManager_ = std::move(other.phaseManager_);
        modelFactory_ = std::move(other.modelFactory_);

        // Rebind raw pointers to new context
        state_ = context_.thermo.get();
        io_ = context_.io.get();
        gemState_ = context_.gem.get();
        phaseConstraints_ = context_.phaseConstraints.get();

        // Rebuild strategies that hold references to state objects
        if (newtonSolver_) {
            newtonSolver_ = std::make_unique<NewtonSolver>(*io_, phaseConstraints_);
        }
        if (phaseManager_) {
            phaseManager_ = std::make_unique<PhaseAssemblageManager>(*state_, *gemState_, *io_);
        }
    }
    return *this;
}

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
    // Pass reference to PhaseConstraints so both solvers use the same instance
    solver_ = std::make_unique<StandardGEMSolver>(*phaseConstraints_);
}

// =========================================================================
// Database Loading
// =========================================================================

void ThermoClass::setThermoFilename(const std::string& filename) {
    io_->cThermoFileName = filename;
}

int ThermoClass::parseCSDataFile() {
    // Use injected parser if available, otherwise use ChemSageParser directly
    if (parser_) {
        int result = parser_->parse(io_->cThermoFileName,
                                   *state_,
                                   *io_,
                                   *context_.parser);
        context_.setInfoThermo(result);
        if (result == 0) {
            context_.allocateFromParser();
        }
    } else {
        // Use ChemSageParser directly
        int result = ChemSageParser::parse(context_);
        context_.setInfoThermo(result);
    }
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
    // Store raw input value and mark for conversion during initialize()
    // This allows setting temperature before or after setting units
    io_->dTemperatureInput = temperature;
    io_->bTemperatureConverted = false;
}

void ThermoClass::setPressure(double pressure) {
    // Store raw input value and mark for conversion during initialize()
    // This allows setting pressure before or after setting units
    io_->dPressureInput = pressure;
    io_->bPressureConverted = false;
}

void ThermoClass::setTemperaturePressure(double temperature, double pressure) {
    // Store raw input values and mark for conversion during initialize()
    // This allows setting values before or after setting units
    io_->dTemperatureInput = temperature;
    io_->dPressureInput = pressure;
    io_->bTemperatureConverted = false;
    io_->bPressureConverted = false;
}

void ThermoClass::setElementMass(int atomicNumber, double mass) {
    if (atomicNumber >= 1 && atomicNumber < Constants::kMaxIsotopes) {
        io_->dElementMass[atomicNumber] = mass;
    }
}

void ThermoClass::setElementMass(const std::string& elementName, double mass) {
    int atomicNum = getAtomicNumber(elementName);
    if (atomicNum > 0) {
        setElementMass(atomicNum, mass);
    }
}

void ThermoClass::setStandardUnits() {
    setUnits("K", "atm", "moles");
}

void ThermoClass::setUnitsSI() {
    // Set SI units: K, Pa, moles
    setUnits("K", "Pa", "moles");
}

void ThermoClass::setUnits(const std::string& tempUnit, const std::string& pressUnit, const std::string& massUnit) {
    setUnitTemperature(tempUnit);
    setUnitPressure(pressUnit);
    setUnitMass(massUnit);
}

void ThermoClass::setUnitTemperature(const std::string& unit) {
    io_->cInputUnitTemperature = unit;
    // Mark temperature for reconversion with new unit
    io_->bTemperatureConverted = false;
}

void ThermoClass::setUnitPressure(const std::string& unit) {
    io_->cInputUnitPressure = unit;
    // Mark pressure for reconversion with new unit
    io_->bPressureConverted = false;
}

void ThermoClass::setUnitMass(const std::string& unit) {
    io_->cInputUnitMass = unit;
}

void ThermoClass::setFuzzyStoich(bool enable) {
    io_->lFuzzyStoich = enable;
}

void ThermoClass::setFuzzyMagnitude(double magnitude) {
    io_->dFuzzMag = magnitude;
}

void ThermoClass::setPrintResultsMode(int mode) {
    io_->iPrintResultsMode = mode;
}

void ThermoClass::setWriteJSON(bool enable) {
    io_->lWriteJSON = enable;
}

void ThermoClass::setHeatCapacityEntropyEnthalpy(bool enable) {
    io_->lHeatCapacityEntropyEnthalpy = enable;
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
    int phaseIdx = state_->getPhaseIndex(phaseName);

    if (phaseIdx < 0 || phaseIdx >= static_cast<int>(phaseConstraints_->solnPhaseConstraints.size())) {
        // Phase not found in solution phases
        return;
    }

    // Clamp target fraction to valid range
    targetFraction = std::max(0.0, std::min(1.0, targetFraction));

    phaseConstraints_->solnPhaseConstraints[phaseIdx].mode = PhaseConstraintMode::Fixed;
    phaseConstraints_->solnPhaseConstraints[phaseIdx].targetFraction = targetFraction;
    phaseConstraints_->solnPhaseConstraints[phaseIdx].lagrangeMultiplier = 0.0;

    // Switch to constrained solver if constraints are active
    if (phaseConstraints_->hasActiveConstraints()) {
        solver_ = std::make_unique<ConstrainedGEMSolver>(*phaseConstraints_);
    }
}

void ThermoClass::setCondPhaseConstraint(const std::string& speciesName, double targetFraction) {
    int speciesIdx = state_->getSpeciesIndex(speciesName);

    if (speciesIdx < 0) {
        return;  // Species not found
    }

    // Pure condensed species are indexed starting from nSpeciesInSolnPhases
    int nSolnSpecies = state_->nSpeciesPhase(state_->nSolnPhasesSys);
    int condIdx = speciesIdx - nSolnSpecies;

    if (condIdx < 0 || condIdx >= static_cast<int>(phaseConstraints_->condPhaseConstraints.size())) {
        return;  // Not a pure condensed species
    }

    // Clamp target fraction to valid range
    targetFraction = std::max(0.0, std::min(1.0, targetFraction));

    phaseConstraints_->condPhaseConstraints[condIdx].mode = PhaseConstraintMode::Fixed;
    phaseConstraints_->condPhaseConstraints[condIdx].targetFraction = targetFraction;
    phaseConstraints_->condPhaseConstraints[condIdx].lagrangeMultiplier = 0.0;

    // Switch to constrained solver if constraints are active
    if (phaseConstraints_->hasActiveConstraints()) {
        solver_ = std::make_unique<ConstrainedGEMSolver>(*phaseConstraints_);
    }
}

void ThermoClass::removePhaseConstraint(const std::string& phaseName) {
    // Try as solution phase first
    int phaseIdx = state_->getPhaseIndex(phaseName);
    if (phaseIdx >= 0 && phaseIdx < static_cast<int>(phaseConstraints_->solnPhaseConstraints.size())) {
        phaseConstraints_->solnPhaseConstraints[phaseIdx].mode = PhaseConstraintMode::None;
        phaseConstraints_->solnPhaseConstraints[phaseIdx].targetFraction = 0.0;
        phaseConstraints_->solnPhaseConstraints[phaseIdx].lagrangeMultiplier = 0.0;

        // Switch back to standard solver if no constraints remain
        if (!phaseConstraints_->hasActiveConstraints()) {
            solver_ = std::make_unique<StandardGEMSolver>(*phaseConstraints_);
        }
        return;
    }

    // Try as pure condensed species
    int speciesIdx = state_->getSpeciesIndex(phaseName);
    if (speciesIdx >= 0) {
        int nSolnSpecies = state_->nSpeciesPhase(state_->nSolnPhasesSys);
        int condIdx = speciesIdx - nSolnSpecies;
        if (condIdx >= 0 && condIdx < static_cast<int>(phaseConstraints_->condPhaseConstraints.size())) {
            phaseConstraints_->condPhaseConstraints[condIdx].mode = PhaseConstraintMode::None;
            phaseConstraints_->condPhaseConstraints[condIdx].targetFraction = 0.0;
            phaseConstraints_->condPhaseConstraints[condIdx].lagrangeMultiplier = 0.0;
        }
    }

    // Switch back to standard solver if no constraints remain
    if (!phaseConstraints_->hasActiveConstraints()) {
        solver_ = std::make_unique<StandardGEMSolver>(*phaseConstraints_);
    }
}

void ThermoClass::clearPhaseConstraints() {
    phaseConstraints_->clear();

    // Switch back to standard solver
    solver_ = std::make_unique<StandardGEMSolver>(*phaseConstraints_);
}

void ThermoClass::setPhaseConstraint(int phaseIndex, bool isSolutionPhase, double targetFraction) {
    // Clamp target fraction to valid range
    targetFraction = std::max(0.0, std::min(1.0, targetFraction));

    if (isSolutionPhase) {
        if (phaseIndex >= 0 && phaseIndex < static_cast<int>(phaseConstraints_->solnPhaseConstraints.size())) {
            phaseConstraints_->solnPhaseConstraints[phaseIndex].mode = PhaseConstraintMode::Fixed;
            phaseConstraints_->solnPhaseConstraints[phaseIndex].targetFraction = targetFraction;
            phaseConstraints_->solnPhaseConstraints[phaseIndex].lagrangeMultiplier = 0.0;
        }
    } else {
        if (phaseIndex >= 0 && phaseIndex < static_cast<int>(phaseConstraints_->condPhaseConstraints.size())) {
            phaseConstraints_->condPhaseConstraints[phaseIndex].mode = PhaseConstraintMode::Fixed;
            phaseConstraints_->condPhaseConstraints[phaseIndex].targetFraction = targetFraction;
            phaseConstraints_->condPhaseConstraints[phaseIndex].lagrangeMultiplier = 0.0;
        }
    }

    // Switch to constrained solver if constraints are active
    if (phaseConstraints_->hasActiveConstraints()) {
        solver_ = std::make_unique<ConstrainedGEMSolver>(*phaseConstraints_);
    }
}

std::pair<double, int> ThermoClass::getPhaseElementFraction(const std::string& phaseName) const {
    // Compute total element moles in system
    double totalElementMoles = 0.0;
    for (int j = 0; j < state_->nElements - state_->nChargedConstraints; ++j) {
        totalElementMoles += state_->dMolesElement(j);
    }

    if (totalElementMoles <= 0.0) {
        return {0.0, -1};  // No elements in system
    }

    // Try as solution phase first
    int phaseIdx = state_->getPhaseIndex(phaseName);
    if (phaseIdx >= 0 && phaseIdx < state_->nSolnPhasesSys) {
        // Check if phase is in the current assemblage
        bool inAssemblage = false;
        for (int i = state_->nElements - state_->nSolnPhases; i < state_->nElements; ++i) {
            if (state_->iAssemblage(i) == -(phaseIdx + 1)) {
                inAssemblage = true;
                break;
            }
        }

        if (!inAssemblage) {
            return {0.0, 0};  // Phase not in assemblage, fraction is zero
        }

        // Sum element moles in this phase
        double phaseElementMoles = 0.0;
        int iFirst = (phaseIdx > 0) ? state_->nSpeciesPhase(phaseIdx) : 0;
        int iLast = state_->nSpeciesPhase(phaseIdx + 1);

        for (int i = iFirst; i < iLast; ++i) {
            double ppm = static_cast<double>(state_->iParticlesPerMole(i));
            for (int j = 0; j < state_->nElements - state_->nChargedConstraints; ++j) {
                phaseElementMoles += state_->dMolesSpecies(i) * state_->dStoichSpecies(i, j) / ppm;
            }
        }

        return {phaseElementMoles / totalElementMoles, 0};
    }

    // Try as pure condensed species
    int speciesIdx = state_->getSpeciesIndex(phaseName);
    if (speciesIdx >= 0) {
        // Check if species is in the current assemblage
        bool inAssemblage = false;
        for (int i = 0; i < state_->nConPhases; ++i) {
            if (state_->iAssemblage(i) == speciesIdx + 1) {
                inAssemblage = true;
                break;
            }
        }

        if (!inAssemblage) {
            return {0.0, 0};  // Species not in assemblage, fraction is zero
        }

        // Sum element moles from this species
        double ppm = static_cast<double>(state_->iParticlesPerMole(speciesIdx));
        double speciesElementMoles = 0.0;
        for (int j = 0; j < state_->nElements - state_->nChargedConstraints; ++j) {
            speciesElementMoles += state_->dMolesSpecies(speciesIdx) * state_->dStoichSpecies(speciesIdx, j) / ppm;
        }

        return {speciesElementMoles / totalElementMoles, 0};
    }

    return {0.0, -1};  // Phase not found
}

bool ThermoClass::arePhaseConstraintsSatisfied() const {
    return phaseConstraints_->areConstraintsSatisfied();
}

void ThermoClass::setConstraintTolerance(double tolerance) {
    phaseConstraints_->constraintTolerance = std::max(1e-10, tolerance);
}

void ThermoClass::setConstraintPenaltyParameter(double rho) {
    double value = std::max(1e-10, rho);
    phaseConstraints_->penaltyParameter = value;
    phaseConstraints_->initialPenaltyParameter = value;
}

void ThermoClass::setConstraintMaxOuterIterations(int maxIter) {
    phaseConstraints_->maxOuterIterations = std::max(1, maxIter);
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
    io_->INFOThermo = 0;

    // Initialize tolerances
    state_->tolerances.initDefaults();

    // Convert input units to internal units (K, atm) only if not already converted
    // This prevents double-conversion on repeated calculate() calls while allowing
    // temperature/pressure to be set before or after unit configuration
    if (!io_->bTemperatureConverted) {
        io_->dTemperature = io_->convertTemperatureToKelvin(io_->dTemperatureInput);
        io_->bTemperatureConverted = true;
    }

    if (!io_->bPressureConverted) {
        io_->dPressure = io_->convertPressureToAtm(io_->dPressureInput);
        io_->bPressureConverted = true;
    }
}

void ThermoClass::checkSystem() {
    // Check temperature (must be positive, not NaN, and within bounds)
    if (std::isnan(io_->dTemperature) || io_->dTemperature <= 0.0 || io_->dTemperature > 10000.0) {
        io_->INFOThermo = ErrorCode::kTemperatureOutOfRange;
        return;
    }

    // Check pressure (must be positive, not NaN, and within bounds)
    if (std::isnan(io_->dPressure) || io_->dPressure <= 0.0 || io_->dPressure > 1e10) {
        io_->INFOThermo = ErrorCode::kPressureOutOfRange;
        return;
    }

    // Check that database is loaded
    if (state_->nSpecies == 0) {
        io_->INFOThermo = ErrorCode::kNoSpeciesInSystem;
        return;
    }

    // Check element masses - must not be NaN or negative
    for (int i = 0; i < Constants::kMaxIsotopes; ++i) {
        if (std::isnan(io_->dElementMass[i]) || io_->dElementMass[i] < 0.0) {
            io_->INFOThermo = ErrorCode::kCompositionOutOfRange;
            return;
        }
    }

    // Check that at least one element has mass
    bool hasElement = false;
    for (int i = 0; i < Constants::kMaxIsotopes; ++i) {
        if (io_->dElementMass[i] > 0.0) {
            hasElement = true;
            break;
        }
    }

    if (!hasElement) {
        io_->INFOThermo = ErrorCode::kCompositionOutOfRange;
        return;
    }
}

void ThermoClass::computeThermoData() {
    auto& parser = *context_.parser;

    double T = io_->dTemperature;
    double R = Constants::kIdealGasConstant;
    double logT = std::log(T);

    // Pre-compute temperature terms for standard Gibbs energy polynomial
    // G = A1 + A2*T + A3*T*ln(T) + A4*T^2 + A5*T^3 + A6/T + additional terms
    double dGibbsCoeff[6];
    dGibbsCoeff[0] = 1.0;              // Multiplies A1 (constant)
    dGibbsCoeff[1] = T;                // Multiplies A2 (linear in T)
    dGibbsCoeff[2] = T * logT;         // Multiplies A3 (T*ln(T))
    dGibbsCoeff[3] = T * T;            // Multiplies A4 (T^2)
    dGibbsCoeff[4] = T * T * T;        // Multiplies A5 (T^3)
    dGibbsCoeff[5] = 1.0 / T;          // Multiplies A6 (1/T)

    // Compute standard Gibbs energy for each species
    for (int i = 0; i < state_->nSpecies; ++i) {
        double G = 0.0;

        // Get number of Gibbs equations for this species
        int nEqs = (i < parser.nGibbsEqSpecies.size()) ? parser.nGibbsEqSpecies(i) : 1;
        if (nEqs <= 0) nEqs = 1;

        // Find the correct temperature range equation
        int iEq = 0;
        for (int eq = 0; eq < nEqs - 1; ++eq) {
            int baseIdx = eq * ParserState::kNumGibbsCoeff;
            double Tmax = parser.dGibbsCoeffSpeciesTemp(i, baseIdx);
            if (T <= Tmax) {
                iEq = eq;
                break;
            }
            iEq = eq + 1;
        }

        // Get base index for this equation's coefficients
        int baseIdx = iEq * ParserState::kNumGibbsCoeff;

        // Compute Gibbs energy from polynomial coefficients
        // Coefficients are stored at indices 1-6 (index 0 is Tmax)
        for (int k = 0; k < 6; ++k) {
            G += parser.dGibbsCoeffSpeciesTemp(i, baseIdx + k + 1) * dGibbsCoeff[k];
        }

        // Add additional terms (up to 3 pairs of coefficient + exponent)
        // Stored at indices 7-8, 9-10, 11-12
        for (int addTerm = 0; addTerm < 3; ++addTerm) {
            int coeffIdx = baseIdx + 7 + addTerm * 2;
            int expIdx = coeffIdx + 1;

            if (expIdx < parser.dGibbsCoeffSpeciesTemp.cols()) {
                double coeff = parser.dGibbsCoeffSpeciesTemp(i, coeffIdx);
                double exp = parser.dGibbsCoeffSpeciesTemp(i, expIdx);

                if (std::abs(coeff) > 1e-30) {
                    if (std::abs(exp - 99.0) < 0.5) {
                        // Exponent 99 means ln(T)
                        G += coeff * logT;
                    } else if (std::abs(exp) > 1e-10) {
                        // General power term: coeff * T^exp
                        G += coeff * std::pow(T, exp);
                    }
                }
            }
        }

        // Store the Gibbs energy
        state_->dStdGibbsEnergy(i) = G;

        // Normalize by RT for chemical potential
        state_->dChemicalPotential(i) = G / (R * T);
    }

    // Compute temperature-dependent excess Gibbs energy parameters
    // G = L0 + L1*T + L2*T*ln(T) + L3*T^2 + L4*T^3 + L5/T
    for (int i = 0; i < state_->nParam; ++i) {
        double L0 = parser.dRegularParamCS(i, 0);
        double L1 = parser.dRegularParamCS(i, 1);
        double L2 = parser.dRegularParamCS(i, 2);
        double L3 = parser.dRegularParamCS(i, 3);
        double L4 = parser.dRegularParamCS(i, 4);
        double L5 = parser.dRegularParamCS(i, 5);

        state_->dExcessGibbsParam(i) = L0 + L1 * T + L2 * T * logT +
                                       L3 * T * T + L4 * T * T * T + L5 / T;
    }

    // Compute element amounts - convert from input mass units to moles
    state_->dMolesElement.setZero();

    // Determine conversion factor based on mass unit (default to "moles" if not set)
    std::string massUnit = io_->cInputUnitMass;
    if (massUnit.empty()) {
        massUnit = "moles";  // Default for backward compatibility
    }
    std::transform(massUnit.begin(), massUnit.end(), massUnit.begin(), ::tolower);

    for (int i = 0; i < Constants::kNumElementsPT; ++i) {
        if (io_->dElementMass[i] > 0.0) {
            // Find element in system
            for (int j = 0; j < state_->nElements; ++j) {
                // Match by atomic number
                int atomicNum = getAtomicNumber(state_->cElementName[j]);
                if (atomicNum == i) {
                    double inputMass = io_->dElementMass[i];
                    double moles = 0.0;

                    // Convert to moles based on input unit
                    if (massUnit == "moles" || massUnit == "mole" || massUnit == "mol") {
                        // Already in moles, no conversion needed
                        moles = inputMass;
                    } else if (massUnit == "grams" || massUnit == "gram" || massUnit == "g") {
                        // Convert grams to moles: moles = grams / atomic_mass
                        double atomicMass = state_->dAtomicMass(j);
                        if (atomicMass <= 0.0) {
                            io_->INFOThermo = 42;  // Invalid atomic mass
                            return;
                        }
                        moles = inputMass / atomicMass;
                    } else if (massUnit == "kilograms" || massUnit == "kilogram" || massUnit == "kg") {
                        // Convert kilograms to moles: moles = kg * 1000 / atomic_mass
                        double atomicMass = state_->dAtomicMass(j);
                        if (atomicMass <= 0.0) {
                            io_->INFOThermo = 42;  // Invalid atomic mass
                            return;
                        }
                        moles = (inputMass * 1000.0) / atomicMass;
                    } else {
                        // Unsupported mass unit
                        io_->INFOThermo = 41;  // Unsupported unit
                        return;
                    }

                    state_->dMolesElement(j) = moles;
                    break;
                }
            }
        }
    }
}

void ThermoClass::setup() {
    // Initialize species moles
    state_->dMolesSpecies.setZero();

    // Run leveling solver to get initial phase assemblage
    levelingSolver(context_);

    if (io_->INFOThermo != 0) return;

    // Run post-leveling solver
    postLevelingSolver(context_);

    // Initialize GEM solver state (already part of solve() - no separate init needed)
}

int ThermoClass::solve() {
    if (!solver_) {
        solver_ = std::make_unique<StandardGEMSolver>(*phaseConstraints_);
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

    ::Thermochimica::postProcess(context_);
    return context_.infoThermo();
}

// =========================================================================
// Output Retrieval
// =========================================================================

std::pair<double, int> ThermoClass::getOutputChemPot(const std::string& elementName) const {
    int elemIdx = state_->getElementIndex(elementName);
    if (elemIdx < 0) {
        return {0.0, -1};
    }

    return {state_->dElementPotential(elemIdx), 0};
}

std::pair<double, int> ThermoClass::getMolesPhase(const std::string& phaseName) const {
    // Try as solution phase first
    int phaseIdx = state_->getPhaseIndex(phaseName);
    if (phaseIdx >= 0) {
        // Find solution phase in assemblage
        for (int i = 0; i < state_->nElements; ++i) {
            if (state_->iAssemblage(i) == -(phaseIdx + 1)) {
                return {state_->dMolesPhase(i), 0};
            }
        }
    }

    // If not found as solution phase, try as pure condensed phase
    int speciesIdx = state_->getSpeciesIndex(phaseName);
    if (speciesIdx >= 0) {
        for (int i = 0; i < state_->nConPhases; ++i) {
            if (state_->iAssemblage(i) == speciesIdx + 1) {
                return {state_->dMolesPhase(i), 0};
            }
        }
    }

    return {0.0, -1};
}

double ThermoClass::getGibbsEnergy() const {
    return io_->dGibbsEnergySys;
}

void ThermoClass::printResults() const {
    // printResults just calls printResultsDetailed
    printResultsDetailed();
}

double ThermoClass::getHeatCapacity() const {
    return io_->dHeatCapacity;
}

double ThermoClass::getEntropy() const {
    return io_->dEntropy;
}

double ThermoClass::getEnthalpy() const {
    return io_->dEnthalpy;
}

std::tuple<double, double, int> ThermoClass::getOutputSolnSpecies(const std::string& phaseName,
                                                                    const std::string& speciesName) const {
    // Get phase index and species range for that phase
    int phaseIdx = state_->getPhaseIndex(phaseName);
    if (phaseIdx < 0) {
        return {0.0, 0.0, -1};  // Phase not found
    }

    // Get species index range for this phase
    int iFirst = (phaseIdx > 0) ? state_->nSpeciesPhase(phaseIdx) : 0;
    int iLast = state_->nSpeciesPhase(phaseIdx + 1);

    // Search for species within this phase's species range
    for (int i = iFirst; i < iLast; ++i) {
        if (state_->cSpeciesName[i] == speciesName) {
            return {state_->dMolFraction(i),
                    state_->dChemicalPotential(i), 0};
        }
    }

    // Species not found in this phase
    return {0.0, 0.0, -1};
}

std::tuple<double, double, int> ThermoClass::getOutputMolSpecies(const std::string& phaseName,
                                                                   const std::string& speciesName) const {
    // Get phase index and species range for that phase
    int phaseIdx = state_->getPhaseIndex(phaseName);
    if (phaseIdx < 0) {
        return {0.0, 0.0, -1};  // Phase not found
    }

    // Get species index range for this phase
    int iFirst = (phaseIdx > 0) ? state_->nSpeciesPhase(phaseIdx) : 0;
    int iLast = state_->nSpeciesPhase(phaseIdx + 1);

    // Search for species within this phase's species range
    for (int i = iFirst; i < iLast; ++i) {
        if (state_->cSpeciesName[i] == speciesName) {
            return {state_->dMolesSpecies(i),
                    state_->dMolFraction(i), 0};
        }
    }

    // Species not found in this phase
    return {0.0, 0.0, -1};
}

std::pair<double, int> ThermoClass::getSolnPhaseMol(const std::string& phaseName) const {
    return getMolesPhase(phaseName);
}

std::pair<double, int> ThermoClass::getPureConPhaseMol(const std::string& phaseName) const {
    // For pure condensed phases, search by species name
    int speciesIdx = state_->getSpeciesIndex(phaseName);
    if (speciesIdx < 0) {
        return {0.0, -1};
    }

    for (int i = 0; i < state_->nConPhases; ++i) {
        if (state_->iAssemblage(i) == speciesIdx + 1) {
            return {state_->dMolesPhase(i), 0};
        }
    }

    return {0.0, -1};
}

std::pair<double, int> ThermoClass::getElementMolesInPhase(const std::string& elementName,
                                                             const std::string& phaseName) const {
    int elemIdx = state_->getElementIndex(elementName);
    int phaseIdx = state_->getPhaseIndex(phaseName);

    if (elemIdx < 0 || phaseIdx < 0) {
        return {0.0, -1};
    }

    // Sum element contribution from all species in phase
    double sum = 0.0;
    int iFirst = (phaseIdx > 0) ? state_->nSpeciesPhase(phaseIdx) : 0;
    int iLast = state_->nSpeciesPhase(phaseIdx + 1);

    for (int i = iFirst; i < iLast; ++i) {
        // Normalize by particles per mole for multi-particle species (e.g., ionic formula units)
        double ppm = static_cast<double>(state_->iParticlesPerMole(i));
        sum += state_->dMolesSpecies(i) * state_->dStoichSpecies(i, elemIdx) / ppm;
    }

    return {sum, 0};
}

int ThermoClass::getPhaseIndex(const std::string& phaseName) const {
    return state_->getPhaseIndex(phaseName);
}

std::pair<double, int> ThermoClass::getOutputSiteFraction(const std::string& phaseName,
                                                           int /*sublattice*/,
                                                           const std::string& /*constituentName*/) const {
    int phaseIdx = state_->getPhaseIndex(phaseName);
    if (phaseIdx < 0) {
        return {0.0, -1};
    }

    int subIdx = state_->iPhaseSublattice(phaseIdx);
    if (subIdx <= 0 || static_cast<size_t>(subIdx - 1) >= state_->dSiteFraction.size()) {
        return {0.0, -1};
    }

    // Find constituent index - simplified (would need to search constituent names)
    return {0.0, -1};
}

bool ThermoClass::isPhaseGas(const std::string& phaseName) const {
    int phaseIdx = state_->getPhaseIndex(phaseName);
    if (phaseIdx < 0 || phaseIdx >= static_cast<int>(state_->iSolnPhaseType.size())) {
        return false;
    }

    Constants::PhaseType type = state_->iSolnPhaseType[phaseIdx];
    // IDMX (ideal mixing) is often used for gas phases
    return type == Constants::PhaseType::IDMX;
}

bool ThermoClass::isPhaseMQM(const std::string& phaseName) const {
    int phaseIdx = state_->getPhaseIndex(phaseName);
    if (phaseIdx < 0 || phaseIdx >= static_cast<int>(state_->iSolnPhaseType.size())) {
        return false;
    }

    Constants::PhaseType type = state_->iSolnPhaseType[phaseIdx];
    return type == Constants::PhaseType::SUBG || type == Constants::PhaseType::SUBQ;
}

std::pair<double, int> ThermoClass::getElementChemicalPotential(int elementIndex) const {
    if (elementIndex < 0 || elementIndex >= state_->nElements) {
        return {0.0, -1};
    }
    return {state_->dElementPotential(elementIndex), 0};
}

std::vector<double> ThermoClass::getAllElementChemicalPotentials() const {
    std::vector<double> potentials(state_->nElements);
    for (int j = 0; j < state_->nElements; ++j) {
        potentials[j] = state_->dElementPotential(j);
    }
    return potentials;
}

std::pair<double, int> ThermoClass::getGibbsEnergyDerivative(const std::string& phaseName) const {
    // Check solution phases first
    int phaseIdx = state_->getPhaseIndex(phaseName);
    if (phaseIdx >= 0 && phaseIdx < static_cast<int>(phaseConstraints_->solnPhaseConstraints.size())) {
        auto& c = phaseConstraints_->solnPhaseConstraints[phaseIdx];
        if (c.mode == PhaseConstraintMode::Fixed) {
            // dG/df = -λ (negative Lagrange multiplier)
            // The Lagrange multiplier is dimensionless (divided by RT internally)
            // To return in J, multiply by RT
            double R = 8.314462618;  // J/(mol·K)
            double T = io_->dTemperature;
            return {-c.lagrangeMultiplier * R * T, 0};
        }
        // Phase exists but no constraint - derivative is 0 at equilibrium
        return {0.0, 0};
    }

    // Check pure condensed phases
    int speciesIdx = state_->getSpeciesIndex(phaseName);
    if (speciesIdx >= 0) {
        int nSolnSpecies = state_->nSpeciesPhase(state_->nSolnPhasesSys);
        int condIdx = speciesIdx - nSolnSpecies;
        if (condIdx >= 0 && condIdx < static_cast<int>(phaseConstraints_->condPhaseConstraints.size())) {
            auto& c = phaseConstraints_->condPhaseConstraints[condIdx];
            if (c.mode == PhaseConstraintMode::Fixed) {
                double R = 8.314462618;
                double T = io_->dTemperature;
                return {-c.lagrangeMultiplier * R * T, 0};
            }
            return {0.0, 0};
        }
    }

    return {0.0, -1};  // Phase not found
}

void ThermoClass::printResultsDetailed() const {
    std::cout << "\n";
    std::cout << "========================================\n";
    std::cout << "        THERMOCHIMICA RESULTS\n";
    std::cout << "========================================\n";

    std::cout << "\nConditions:\n";
    std::cout << "  Temperature: " << std::fixed << std::setprecision(2)
              << io_->dTemperature << " K\n";
    std::cout << "  Pressure:    " << io_->dPressure << " atm\n";

    std::cout << "\nSystem Gibbs Energy: " << std::scientific << std::setprecision(6)
              << io_->dGibbsEnergySys << " J\n";

    if (state_->nSolnPhases > 0) {
        std::cout << "\nStable Solution Phases:\n";
        for (int i = 0; i < state_->nSolnPhases; ++i) {
            int idx = state_->nElements - state_->nSolnPhases + i;
            int phaseIdx = -state_->iAssemblage(idx) - 1;
            if (phaseIdx >= 0) {
                std::cout << "  " << std::setw(20) << std::left
                          << state_->cSolnPhaseName[phaseIdx]
                          << ": " << std::scientific << std::setprecision(4)
                          << state_->dMolesPhase(idx) << " mol\n";
            }
        }
    }

    if (state_->nConPhases > 0) {
        std::cout << "\nStable Pure Condensed Phases:\n";
        for (int i = 0; i < state_->nConPhases; ++i) {
            int speciesIdx = state_->iAssemblage(i) - 1;
            if (speciesIdx >= 0) {
                std::cout << "  " << std::setw(20) << std::left
                          << state_->cSpeciesName[speciesIdx]
                          << ": " << std::scientific << std::setprecision(4)
                          << state_->dMolesPhase(i) << " mol\n";
            }
        }
    }

    std::cout << "\n========================================\n\n";
}

void ThermoClass::writeJSON(bool append) {
    (void)append;
    static bool warned = false;
    if (!warned) {
        std::cerr << "Warning: writeJSON is not yet implemented\n";
        warned = true;
    }
}

void ThermoClass::postProcess() {
    ::Thermochimica::postProcess(context_);
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
    return ErrorCode::getMessage(context_.infoThermo());
}

// =========================================================================
// Database Queries
// =========================================================================

int ThermoClass::getNumberElementsDatabase() const {
    return state_->nElements;
}

std::string ThermoClass::getElementAtIndex(int index) const {
    if (index >= 0 && index < static_cast<int>(state_->cElementName.size())) {
        return state_->cElementName[index];
    }
    return "";
}

std::pair<int, int> ThermoClass::getNumberPhasesSystem() const {
    return {state_->nSolnPhasesSys, state_->nConPhasesSys};
}

std::string ThermoClass::getPhaseNameAtIndex(int index) const {
    if (index >= 0 && index < static_cast<int>(state_->cSolnPhaseName.size())) {
        return state_->cSolnPhaseName[index];
    }
    return "";
}

int ThermoClass::getNumberSpeciesSystem() const {
    return state_->nSpecies;
}

std::string ThermoClass::getSpeciesAtIndex(int index) const {
    if (index >= 0 && index < static_cast<int>(state_->cSpeciesName.size())) {
        return state_->cSpeciesName[index];
    }
    return "";
}

// =========================================================================
// Reinitialization
// =========================================================================

void ThermoClass::saveReinitData() {
    auto& reinit = *context_.reinit;

    reinit.saveState(state_->iAssemblage,
                     state_->dChemicalPotential,
                     state_->dMolesPhase,
                     state_->dElementPotential,
                     state_->dMolFraction);

    io_->lReinitAvailable = true;
}

void ThermoClass::setReinitRequested(bool requested) {
    io_->lReinitRequested = requested;
}

bool ThermoClass::isReinitDataAvailable() const {
    return io_->lReinitAvailable && context_.reinit->isValid();
}

std::tuple<std::vector<int>, std::vector<double>, std::vector<double>,
           std::vector<double>, std::vector<double>> ThermoClass::getReinitData() const {
    auto& reinit = *context_.reinit;

    std::vector<int> assemblage(reinit.iAssemblage_Old.data(),
                                reinit.iAssemblage_Old.data() + reinit.iAssemblage_Old.size());
    std::vector<double> chemPot(reinit.dChemicalPotential_Old.data(),
                                reinit.dChemicalPotential_Old.data() + reinit.dChemicalPotential_Old.size());
    std::vector<double> molesPhase(reinit.dMolesPhase_Old.data(),
                                   reinit.dMolesPhase_Old.data() + reinit.dMolesPhase_Old.size());
    std::vector<double> elemPot(reinit.dElementPotential_Old.data(),
                                reinit.dElementPotential_Old.data() + reinit.dElementPotential_Old.size());
    std::vector<double> molFrac(reinit.dMolFraction_Old.data(),
                                reinit.dMolFraction_Old.data() + reinit.dMolFraction_Old.size());

    return {assemblage, chemPot, molesPhase, elemPot, molFrac};
}

void ThermoClass::setReinitData(const std::vector<int>& assemblage,
                                  const std::vector<double>& chemPot,
                                  const std::vector<double>& molesPhase,
                                  const std::vector<double>& elemPot,
                                  const std::vector<double>& molFrac) {
    auto& reinit = *context_.reinit;

    reinit.iAssemblage_Old = Eigen::Map<const Eigen::VectorXi>(assemblage.data(), assemblage.size());
    reinit.dChemicalPotential_Old = Eigen::Map<const Eigen::VectorXd>(chemPot.data(), chemPot.size());
    reinit.dMolesPhase_Old = Eigen::Map<const Eigen::VectorXd>(molesPhase.data(), molesPhase.size());
    reinit.dElementPotential_Old = Eigen::Map<const Eigen::VectorXd>(elemPot.data(), elemPot.size());
    reinit.dMolFraction_Old = Eigen::Map<const Eigen::VectorXd>(molFrac.data(), molFrac.size());

    io_->lReinitAvailable = true;
}

// =========================================================================
// Reset
// =========================================================================

void ThermoClass::reset() {
    context_.resetThermo();
    // Mark values for reconversion after reset
    io_->bTemperatureConverted = false;
    io_->bPressureConverted = false;
}

void ThermoClass::resetAll() {
    context_.resetAll();
    initializeDefaultStrategies();
    // Mark values for reconversion after full reset
    io_->bTemperatureConverted = false;
    io_->bPressureConverted = false;
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

// =========================================================================
// Utility Functions (Static)
// =========================================================================

int ThermoClass::getAtomicNumber(const std::string& symbol) {
    // Element symbols for lookup
    static const char* kElementSymbols[] = {
        "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
        "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
        "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
    };

    std::string upper = symbol;
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);

    // Handle special cases
    if (upper == "E-" || upper == "E(-)" || upper.substr(0, 2) == "E(") {
        return 0;  // Electron
    }

    // Search element symbols
    for (int i = 1; i <= 118; ++i) {
        std::string elemUpper = kElementSymbols[i];
        std::transform(elemUpper.begin(), elemUpper.end(), elemUpper.begin(), ::toupper);
        if (upper == elemUpper) {
            return i;
        }
    }

    return -1;  // Not found
}

std::string ThermoClass::getElementSymbol(int atomicNumber) {
    static const char* kElementSymbols[] = {
        "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
        "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
        "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
    };

    if (atomicNumber >= 1 && atomicNumber <= 118) {
        return kElementSymbols[atomicNumber];
    }
    return "";
}

} // namespace Thermochimica
