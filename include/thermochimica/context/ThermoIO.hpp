#pragma once

#include <Eigen/Dense>
#include <string>
#include <array>
#include <vector>
#include "../util/Constants.hpp"

namespace Thermochimica {

/// Input/Output state - converted from ModuleThermoIO.f90
/// Bridge between external callers and internal state
struct ThermoIO {
    // Input variables
    int iCounter = 0;                   ///< Counter variable
    int iPrintResultsMode = 0;          ///< Print results mode (0=none, 1=basic, 2=detailed)
    int nMinSpeciesPerPhase = 2;        ///< Minimum species threshold per phase

    double dTemperature = 0.0;          ///< Temperature (internal: K, converted)
    double dPressure = 0.0;             ///< Pressure (internal: atm, converted)
    double dTemperatureInput = 0.0;     ///< Temperature (raw user input, in user units)
    double dPressureInput = 0.0;        ///< Pressure (raw user input, in user units)
    bool bTemperatureConverted = false; ///< Temperature converted for current cycle
    bool bPressureConverted = false;    ///< Pressure converted for current cycle
    double dFuzzMag = 1.0e-12;          ///< Fuzzy stoichiometry magnitude

    std::array<double, Constants::kMaxIsotopes> dElementMass{};  ///< Element masses (0:168)
    std::array<bool, Constants::kNumElementsPT + 1> lPreset{};   ///< Preset element flags (0:118)

    std::string cInputUnitTemperature;  ///< Temperature unit ('K', 'C', 'F', 'R')
    std::string cInputUnitPressure;     ///< Pressure unit ('atm', 'psi', 'bar', 'Pa', 'kPa')
    std::string cInputUnitMass;         ///< Mass unit ('moles', 'grams', etc.)
    std::string cThermoFileName;        ///< Thermodynamic data file path
    std::string cOutputFilePath;        ///< Output file path
    std::string cResolvedOutputFilePath; ///< Resolved output path

    // Control flags
    bool lReinitAvailable = false;      ///< Reinit data available
    bool lReinitLoaded = false;         ///< Reinit data loaded
    bool lReinitRequested = false;      ///< Reinit requested
    bool lStepTogether = false;         ///< Step together flag
    bool lWriteJSON = false;            ///< Write JSON output
    bool lFuzzyStoich = false;          ///< Fuzzy stoichiometry enabled
    bool lGibbsMinCheck = false;        ///< Gibbs minimum check

    // Phase exclusion
    int nPhasesExcluded = 0;            ///< Number of excluded phases
    int nPhasesExcludedExcept = 0;      ///< Number of except phases
    std::vector<std::string> cPhasesExcluded;       ///< Excluded phase names
    std::vector<std::string> cPhasesExcludedExcept; ///< Except phase names

    // Compound variables
    int nCompounds = 0;                 ///< Number of compounds
    std::array<double, Constants::kNumElementsPT> dCompoundMass{};  ///< Compound masses
    Eigen::Matrix<double, Constants::kNumElementsPT, Constants::kNumElementsPT + 1> dCompoundStoich;
    std::array<std::string, Constants::kNumElementsPT> cCompoundNames;
    bool lCompoundStoichCalculated = false;
    bool lRetryAttempted = false;
    bool lHeatCapacityEntropyEnthalpy = false;

    // Output variables
    int INFOThermo = 0;                 ///< Status/error code
    int nSolnPhasesOut = 0;             ///< Output solution phases count
    int nPureConPhaseOut = 0;           ///< Output pure condensed phases count
    int nSpeciesOut = 0;                ///< Output species count
    double dGibbsEnergySys = 0.0;       ///< System Gibbs energy
    double dHeatCapacity = 0.0;         ///< Heat capacity
    double dEntropy = 0.0;              ///< Entropy
    double dEnthalpy = 0.0;             ///< Enthalpy

    std::vector<double> dSolnPhaseMolesOut;       ///< Solution phase moles output
    std::vector<double> dPureConPhaseMolesOut;   ///< Pure condensed phase moles output
    std::vector<double> dSpeciesMoleFractionOut; ///< Species mole fractions output
    std::vector<std::string> cSolnPhaseNameOut;   ///< Solution phase names output
    std::vector<std::string> cPureConPhaseNameOut; ///< Pure condensed phase names output
    std::vector<std::string> cSpeciesNameOut;     ///< Species names output
    std::vector<std::string> cSpeciesPhaseOut;    ///< Species phase assignment output
    std::vector<bool> lSpeciesStable;             ///< Species stability flags

    // Constructor
    ThermoIO();

    /// Check if path is absolute
    static bool pathIsAbsolute(const std::string& path);

    /// Set default output file path
    void setDefaultOutputFilePath();

    /// Update output file path
    void updateOutputFilePath(const std::string& rawPath);

    /// Get resolved output file path
    std::string getResolvedOutputFilePath() const;

    /// Convert temperature to Kelvin
    double convertTemperatureToKelvin(double temp) const;

    /// Convert pressure to atm
    double convertPressureToAtm(double pressure) const;

    /// Reset I/O state
    void reset();

    /// Reset all element masses to zero
    void resetElementMasses();
};

} // namespace Thermochimica
