#pragma once

namespace Thermochimica {
namespace ErrorCode {

// Success
constexpr int kSuccess = 0;

// Input validation errors (1-9)
constexpr int kTemperatureOutOfRange = 1;
constexpr int kPressureOutOfRange = 2;
constexpr int kCompositionOutOfRange = 3;
constexpr int kUnrecognizedTemperatureUnit = 4;
constexpr int kUnrecognizedPressureUnit = 5;
constexpr int kNoDataFileSpecified = 6;
constexpr int kUnrecognizedMassUnit = 7;
constexpr int kInvalidCompoundStoichiometry = 8;

// Parser errors (10-19)
constexpr int kDataFileNotFound = 10;
constexpr int kParserError = 11;
constexpr int kUnsupportedPhaseType = 12;
constexpr int kInvalidDataFormat = 13;

// System check errors (20-29)
constexpr int kNoElementsInSystem = 20;
constexpr int kNoSpeciesInSystem = 21;
constexpr int kSystemCheckFailed = 22;
constexpr int kInsufficientElements = 23;

// Solver errors (30-49)
constexpr int kGEMSolverDidNotConverge = 30;
constexpr int kLevelingSolverFailed = 31;
constexpr int kSingularMatrix = 32;
constexpr int kNoStablePhases = 33;
constexpr int kMaxIterationsReached = 34;
constexpr int kNumericalInstability = 35;
constexpr int kPhaseRuleViolation = 36;

// Post-processing errors (50-59)
constexpr int kPostProcessError = 50;
constexpr int kOutputError = 51;
constexpr int kJSONWriteError = 52;

// Reinitialization errors (60-69)
constexpr int kReinitDataNotAvailable = 60;
constexpr int kReinitDataInvalid = 61;

// Get error message string
inline const char* getMessage(int code) {
    switch (code) {
        case kSuccess: return "Success";
        case kTemperatureOutOfRange: return "Temperature out of range";
        case kPressureOutOfRange: return "Pressure out of range";
        case kCompositionOutOfRange: return "Composition out of range";
        case kUnrecognizedTemperatureUnit: return "Unrecognized temperature unit";
        case kUnrecognizedPressureUnit: return "Unrecognized pressure unit";
        case kNoDataFileSpecified: return "No thermodynamic data file specified";
        case kUnrecognizedMassUnit: return "Unrecognized mass unit";
        case kInvalidCompoundStoichiometry: return "Invalid compound stoichiometry";
        case kDataFileNotFound: return "Data file not found";
        case kParserError: return "Parser error";
        case kUnsupportedPhaseType: return "Unsupported phase type";
        case kInvalidDataFormat: return "Invalid data format";
        case kNoElementsInSystem: return "No elements in system";
        case kNoSpeciesInSystem: return "No species in system";
        case kSystemCheckFailed: return "System check failed";
        case kInsufficientElements: return "Insufficient elements for calculation";
        case kGEMSolverDidNotConverge: return "GEM solver did not converge";
        case kLevelingSolverFailed: return "Leveling solver failed";
        case kSingularMatrix: return "Singular matrix encountered";
        case kNoStablePhases: return "No stable phases found";
        case kMaxIterationsReached: return "Maximum iterations reached";
        case kNumericalInstability: return "Numerical instability detected";
        case kPhaseRuleViolation: return "Phase rule violation";
        case kPostProcessError: return "Post-processing error";
        case kOutputError: return "Output error";
        case kJSONWriteError: return "JSON write error";
        case kReinitDataNotAvailable: return "Reinitialization data not available";
        case kReinitDataInvalid: return "Reinitialization data invalid";
        default: return "Unknown error";
    }
}

} // namespace ErrorCode
} // namespace Thermochimica
