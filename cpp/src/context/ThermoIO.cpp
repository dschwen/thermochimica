#include "thermochimica/context/ThermoIO.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>

namespace Thermochimica {

ThermoIO::ThermoIO() {
    dElementMass.fill(0.0);
    lPreset.fill(false);
    dCompoundMass.fill(0.0);
    dCompoundStoich.setZero();

    cPhasesExcluded.reserve(100);
    cPhasesExcludedExcept.reserve(100);

    setDefaultOutputFilePath();
}

bool ThermoIO::pathIsAbsolute(const std::string& path) {
    if (path.empty()) {
        return false;
    }

    // Unix absolute path
    if (path[0] == '/' || path[0] == '\\') {
        return true;
    }

    // Relative path starting with ..
    if (path.size() >= 2 && path.substr(0, 2) == "..") {
        return true;
    }

    // Windows drive letter (e.g., C:)
    if (path.size() >= 2 && std::isalpha(path[0]) && path[1] == ':') {
        return true;
    }

    return false;
}

void ThermoIO::setDefaultOutputFilePath() {
    updateOutputFilePath("../outputs/thermoout.json");
}

void ThermoIO::updateOutputFilePath(const std::string& rawPath) {
    std::string cleaned = rawPath;

    // Trim whitespace
    auto start = cleaned.find_first_not_of(" \t\n\r");
    auto end = cleaned.find_last_not_of(" \t\n\r");

    if (start == std::string::npos) {
        setDefaultOutputFilePath();
        return;
    }

    cleaned = cleaned.substr(start, end - start + 1);
    cOutputFilePath = cleaned;

    if (pathIsAbsolute(cleaned)) {
        cResolvedOutputFilePath = cleaned;
    } else if (cleaned.find('/') == std::string::npos &&
               cleaned.find('\\') == std::string::npos) {
        // Just a filename, put in outputs directory
        #ifdef DATA_DIRECTORY
        cResolvedOutputFilePath = std::string(DATA_DIRECTORY) + "../outputs/" + cleaned;
        #else
        cResolvedOutputFilePath = "../outputs/" + cleaned;
        #endif
    } else {
        #ifdef DATA_DIRECTORY
        cResolvedOutputFilePath = std::string(DATA_DIRECTORY) + cleaned;
        #else
        cResolvedOutputFilePath = cleaned;
        #endif
    }
}

std::string ThermoIO::getResolvedOutputFilePath() const {
    if (cResolvedOutputFilePath.empty()) {
        return "../outputs/thermoout.json";
    }
    return cResolvedOutputFilePath;
}

double ThermoIO::convertTemperatureToKelvin(double temp) const {
    std::string unit = cInputUnitTemperature;
    std::transform(unit.begin(), unit.end(), unit.begin(), ::toupper);

    if (unit == "K" || unit == "KELVIN" || unit.empty()) {
        return temp;
    } else if (unit == "C" || unit == "CELSIUS") {
        return temp + 273.15;
    } else if (unit == "F" || unit == "FAHRENHEIT") {
        return (temp - 32.0) * 5.0 / 9.0 + 273.15;
    } else if (unit == "R" || unit == "RANKINE") {
        return temp * 5.0 / 9.0;
    }

    // Unknown unit, assume Kelvin
    return temp;
}

double ThermoIO::convertPressureToAtm(double pressure) const {
    std::string unit = cInputUnitPressure;
    std::transform(unit.begin(), unit.end(), unit.begin(), ::toupper);

    if (unit == "ATM" || unit == "ATMOSPHERES" || unit.empty()) {
        return pressure;
    } else if (unit == "PA" || unit == "PASCAL" || unit == "PASCALS") {
        return pressure / 101325.0;
    } else if (unit == "KPA" || unit == "KILOPASCAL" || unit == "KILOPASCALS") {
        return pressure / 101.325;
    } else if (unit == "BAR" || unit == "BARS") {
        return pressure / 1.01325;
    } else if (unit == "PSI") {
        return pressure / 14.6959;
    }

    // Unknown unit, assume atm
    return pressure;
}

void ThermoIO::reset() {
    iCounter = 0;
    iPrintResultsMode = 0;
    nMinSpeciesPerPhase = 2;

    dTemperature = 0.0;
    dPressure = 0.0;
    dFuzzMag = 1.0e-12;

    resetElementMasses();
    lPreset.fill(false);

    cInputUnitTemperature.clear();
    cInputUnitPressure.clear();
    cInputUnitMass.clear();
    cThermoFileName.clear();

    lReinitAvailable = false;
    lReinitLoaded = false;
    lReinitRequested = false;
    lStepTogether = false;
    lWriteJSON = false;
    lFuzzyStoich = false;
    lGibbsMinCheck = false;

    nPhasesExcluded = 0;
    nPhasesExcludedExcept = 0;
    cPhasesExcluded.clear();
    cPhasesExcludedExcept.clear();

    nCompounds = 0;
    dCompoundMass.fill(0.0);
    dCompoundStoich.setZero();
    lCompoundStoichCalculated = false;
    lRetryAttempted = false;
    lHeatCapacityEntropyEnthalpy = false;

    INFOThermo = 0;
    nSolnPhasesOut = 0;
    nPureConPhaseOut = 0;
    nSpeciesOut = 0;
    dGibbsEnergySys = 0.0;
    dHeatCapacity = 0.0;
    dEntropy = 0.0;
    dEnthalpy = 0.0;

    dSolnPhaseMolesOut.clear();
    dPureConPhaseMolesOut.clear();
    dSpeciesMoleFractionOut.clear();
    cSolnPhaseNameOut.clear();
    cPureConPhaseNameOut.clear();
    cSpeciesNameOut.clear();
    cSpeciesPhaseOut.clear();
    lSpeciesStable.clear();

    setDefaultOutputFilePath();
}

void ThermoIO::resetElementMasses() {
    dElementMass.fill(0.0);
}

} // namespace Thermochimica
