#pragma once

namespace Thermochimica {
namespace Constants {

// Physical constants
constexpr double kIdealGasConstant = 8.3144598;  // J/(mol·K)
constexpr double kBoltzmann = 1.38064852e-23;     // J/K
constexpr double kAvogadro = 6.022140857e23;      // 1/mol

// System limits
constexpr int kNumElementsPT = 118;               // Periodic table elements
constexpr int kMaxIsotopes = 169;                 // Including isotope variants (0:168)
constexpr int kMaxCompounds = 50;                 // Maximum compounds
constexpr int kNumTolerances = 15;                // Number of tolerance values
constexpr int kMaxSublattices = 5;                // Maximum sublattices per phase
constexpr int kMaxSolnPhases = 100;               // Maximum solution phases
constexpr int kIterGlobalMax = 3000;              // Maximum GEM iterations

// Gibbs energy coefficients
constexpr int kNumGibbsCoeff = 13;                // Gibbs polynomial coefficients
constexpr int kMaxGibbsEqs = 10;                  // Max Gibbs equations per species
constexpr int kMaxParams = 4;                     // Max mixing parameters

// Numerical constants
constexpr double kMinMoleFraction = 1.0e-100;     // Minimum mole fraction
constexpr double kDefaultNormalize = 1000.0;      // Default normalization factor

// Supported solution phase types
enum class PhaseType {
    Unknown = 0,
    IDMX,       // Ideal mixing
    QKTO,       // Kohler-Toop
    SUBL,       // Compound Energy Formalism (sublattice)
    RKMP,       // Redlich-Kister-Muggianu Polynomial
    RKMPM,      // RKMP with magnetic
    SUBLM,      // SUBL with magnetic
    SUBG,       // Modified Quasichemical Model (MQM) - quadruplets
    SUBQ,       // MQM variant
    SUBI,       // Ionic sublattice
    SUBM        // Sublattice with miscibility
};

// Phase type names for parsing
constexpr const char* kPhaseTypeNames[] = {
    "UNKNOWN",
    "IDMX",
    "QKTO",
    "SUBL",
    "RKMP",
    "RKMPM",
    "SUBLM",
    "SUBG",
    "SUBQ",
    "SUBI",
    "SUBM"
};

constexpr int kNumPhaseTypes = 11;

} // namespace Constants
} // namespace Thermochimica
