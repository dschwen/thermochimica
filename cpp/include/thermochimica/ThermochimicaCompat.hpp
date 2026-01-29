/// @file ThermochimicaCompat.hpp
/// @brief Compatibility layer for migrating from Fortran-backed C++ API
/// @details Provides a global context and wrapper functions that match the
/// old API signatures. This allows gradual migration with minimal code changes.
///
/// Usage:
///   #include <thermochimica/ThermochimicaCompat.hpp>
///   using namespace Thermochimica::Compat;
///
///   // Old-style API calls now work:
///   setThermoFilename("database.dat");
///   parseThermoFile();
///   thermochimica();
///
/// WARNING: This compatibility layer uses a global context and is NOT thread-safe.
/// For thread-safe operation, use the new context-based API directly.

#pragma once

#include "Thermochimica.hpp"
#include "ThermoContext.hpp"

namespace Thermochimica {
namespace Compat {

/// Get the global context (for advanced use)
/// @warning Not thread-safe
inline ThermoContext& globalContext() {
    static ThermoContext ctx;
    return ctx;
}

// ============================================================================
// Setup Functions
// ============================================================================

inline void setThermoFilename(const std::string& filename) {
    Thermochimica::setThermoFilename(globalContext(), filename);
}

inline void parseThermoFile() {
    Thermochimica::parseCSDataFile(globalContext());
}

inline void setStandardUnits() {
    Thermochimica::setStandardUnits(globalContext());
}

inline void setModelicaUnits() {
    Thermochimica::setUnits(globalContext(), "K", "Pa", "kg/m3");
}

inline void setUnitTemperature(const std::string& unit) {
    Thermochimica::setUnitTemperature(globalContext(), unit);
}

inline void setUnitPressure(const std::string& unit) {
    Thermochimica::setUnitPressure(globalContext(), unit);
}

inline void setUnitMass(const std::string& unit) {
    Thermochimica::setUnitMass(globalContext(), unit);
}

inline void setUnits(const std::string& tunit,
                     const std::string& punit,
                     const std::string& munit) {
    Thermochimica::setUnits(globalContext(), tunit, punit, munit);
}

inline void setTemperaturePressure(double temperature, double pressure) {
    Thermochimica::setTemperaturePressure(globalContext(), temperature, pressure);
}

inline void setElementMass(int atomicNumber, double mass) {
    Thermochimica::setElementMass(globalContext(), atomicNumber, mass);
}

inline void presetElementMass(int atomicNumber, double mass) {
    Thermochimica::setElementMass(globalContext(), atomicNumber, mass);
}

// ============================================================================
// Solver Functions
// ============================================================================

inline void thermochimica() {
    Thermochimica::thermochimica(globalContext());
}

inline void setup() {
    Thermochimica::setup(globalContext());
}

inline void solve() {
    Thermochimica::solve(globalContext());
}

inline void init() {
    Thermochimica::init(globalContext());
}

inline void checkSystem() {
    Thermochimica::checkSystem(globalContext());
}

inline void compThermoData() {
    Thermochimica::compThermoData(globalContext());
}

// ============================================================================
// Status and Reset Functions
// ============================================================================

inline int checkInfoThermo() {
    return globalContext().infoThermo();
}

inline void resetInfoThermo() {
    globalContext().io->INFOThermo = 0;
}

inline void resetThermo() {
    Thermochimica::resetThermo(globalContext());
}

inline void resetThermoAll() {
    Thermochimica::resetThermoAll(globalContext());
}

// ============================================================================
// Output Functions
// ============================================================================

inline void printResults() {
    Thermochimica::printResults(globalContext());
}

inline void setPrintResultsMode(int mode) {
    Thermochimica::setPrintResultsMode(globalContext(), mode);
}

// ============================================================================
// Data Retrieval Functions
// ============================================================================

inline std::pair<double, int> getOutputChemPot(const std::string& elementName) {
    return Thermochimica::getOutputChemPot(globalContext(), elementName);
}

inline std::tuple<double, double, int>
getOutputSolnSpecies(const std::string& phaseName, const std::string& speciesName) {
    return Thermochimica::getOutputSolnSpecies(globalContext(), phaseName, speciesName);
}

inline std::tuple<double, double, int>
getOutputMolSpecies(const std::string& speciesName) {
    // Old API didn't require phase name, search all phases
    return Thermochimica::getOutputMolSpecies(globalContext(), "", speciesName);
}

inline std::pair<double, int> getSolnPhaseMol(const std::string& phaseName) {
    return Thermochimica::getSolnPhaseMol(globalContext(), phaseName);
}

inline std::pair<double, int> getPureConPhaseMol(const std::string& phaseName) {
    return Thermochimica::getPureConPhaseMol(globalContext(), phaseName);
}

inline std::pair<double, int>
getElementMolesInPhase(const std::string& elementName, const std::string& phaseName) {
    return Thermochimica::getElementMolesInPhase(globalContext(), elementName, phaseName);
}

inline int getPhaseIndex(const std::string& phaseName) {
    return Thermochimica::getPhaseIndex(globalContext(), phaseName);
}

inline std::pair<double, int>
getOutputSiteFraction(const std::string& phaseName, int sublattice, int constituent) {
    // Old API used constituent index, new API uses name
    // This is a simplified compatibility - may need adjustment
    return Thermochimica::getOutputSiteFraction(globalContext(), phaseName, sublattice,
                                                 std::to_string(constituent));
}

// ============================================================================
// Phase Type Query Functions
// ============================================================================

inline bool isPhaseGas(int phaseIndex) {
    std::string phaseName = Thermochimica::getPhaseNameAtIndex(globalContext(), phaseIndex);
    return Thermochimica::isPhaseGas(globalContext(), phaseName);
}

inline bool isPhaseMQM(int phaseIndex) {
    std::string phaseName = Thermochimica::getPhaseNameAtIndex(globalContext(), phaseIndex);
    return Thermochimica::isPhaseMQM(globalContext(), phaseName);
}

// ============================================================================
// Database Query Functions
// ============================================================================

inline std::size_t getNumberElementsDatabase() {
    return static_cast<std::size_t>(
        Thermochimica::getNumberElementsDatabase(globalContext()));
}

inline std::string getElementAtIndex(int index) {
    return Thermochimica::getElementAtIndex(globalContext(), index);
}

inline std::vector<std::string> getElementsDatabase() {
    auto n = getNumberElementsDatabase();
    std::vector<std::string> elements(n);
    for (std::size_t i = 0; i < n; ++i) {
        elements[i] = getElementAtIndex(static_cast<int>(i));
    }
    return elements;
}

inline std::pair<std::size_t, std::size_t> getNumberPhasesSystem() {
    auto [soln, cond] = Thermochimica::getNumberPhasesSystem(globalContext());
    return {static_cast<std::size_t>(soln), static_cast<std::size_t>(cond)};
}

inline std::string getPhaseNameAtIndex(int index) {
    return Thermochimica::getPhaseNameAtIndex(globalContext(), index);
}

inline std::vector<std::string> getPhaseNamesSystem() {
    auto [nSoln, nCond] = getNumberPhasesSystem();
    auto n = nSoln + nCond;
    std::vector<std::string> names(n);
    for (std::size_t i = 0; i < n; ++i) {
        names[i] = getPhaseNameAtIndex(static_cast<int>(i));
    }
    return names;
}

inline int getNumberSpeciesSystem() {
    return Thermochimica::getNumberSpeciesSystem(globalContext());
}

inline std::string getSpeciesAtIndex(int index) {
    return Thermochimica::getSpeciesAtIndex(globalContext(), index);
}

// ============================================================================
// Thermodynamic Properties
// ============================================================================

inline double getGibbsEnergy() {
    return Thermochimica::getGibbsEnergy(globalContext());
}

inline double getHeatCapacity() {
    return Thermochimica::getHeatCapacity(globalContext());
}

inline double getEntropy() {
    return Thermochimica::getEntropy(globalContext());
}

inline double getEnthalpy() {
    return Thermochimica::getEnthalpy(globalContext());
}

inline void setHeatCapacityEnthalpyEntropyRequested(bool requested) {
    Thermochimica::setHeatCapacityEntropyEnthalpy(globalContext(), requested);
}

inline std::tuple<double, double, double> getHeatCapacityEnthalpyEntropy() {
    return {getHeatCapacity(), getEnthalpy(), getEntropy()};
}

// ============================================================================
// Reinitialization Functions
// ============================================================================

inline void saveReinitData() {
    Thermochimica::saveReinitData(globalContext());
}

inline void setReinitRequested(bool requested) {
    Thermochimica::setReinitRequested(globalContext(), requested);
}

inline bool isReinitDataAvailable() {
    return Thermochimica::isReinitDataAvailable(globalContext());
}

inline std::pair<int, int> getReinitDataSizes() {
    auto& ctx = globalContext();
    return {ctx.thermo->nElements, ctx.thermo->nSpecies};
}

inline std::vector<double> getMolesPhase() {
    auto& thermo = *globalContext().thermo;
    return std::vector<double>(thermo.dMolesPhase.data(),
                               thermo.dMolesPhase.data() + thermo.nElements);
}

inline std::vector<int> getAssemblage() {
    auto& thermo = *globalContext().thermo;
    return std::vector<int>(thermo.iAssemblage.data(),
                            thermo.iAssemblage.data() + thermo.nElements);
}

inline std::vector<double> getAllElementPotential() {
    auto& thermo = *globalContext().thermo;
    return std::vector<double>(thermo.dElementPotential.data(),
                               thermo.dElementPotential.data() + thermo.nElements);
}

// ============================================================================
// Advanced Functions
// ============================================================================

inline void setFuzzyStoich(bool enable) {
    Thermochimica::setFuzzyStoich(globalContext(), enable);
}

inline void setFuzzyMagnitude(double magnitude) {
    Thermochimica::setFuzzyMagnitude(globalContext(), magnitude);
}

inline void setGibbsMinCheck(bool /*requested*/) {
    // Always performed in new API - no-op for compatibility
}

}} // namespace Thermochimica::Compat
