#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <array>
#include "../util/Constants.hpp"

namespace Thermochimica {

/// Parser temporary storage - converted from ModuleParseCS.f90
/// Temporary data structures for ChemSage file parsing
struct ParserState {
    // Constants
    static constexpr int kSolnPhasesSysMax = 100;
    static constexpr int kMaxSublattice = 5;
    static constexpr int kNumSolnTypeSupport = 10;
    static constexpr int kNumGibbsCoeff = 13;
    static constexpr int kMaxGibbsEqs = 6;
    static constexpr int kParamMax = 4;

    // Supported phase type names
    static constexpr const char* kSolnPhaseTypeSupport[kNumSolnTypeSupport] = {
        "IDMX", "QKTO", "SUBL", "RKMP", "RKMPM",
        "SUBLM", "SUBG", "SUBQ", "SUBI", "SUBM"
    };

    // Scalar integers
    int nElementsCS = 0;            ///< Number of elements parsed
    int nSpeciesCS = 0;             ///< Number of species parsed
    int nSolnPhasesSysCS = 0;       ///< Number of solution phases
    int INFO = 0;                   ///< Parser status
    int iMiscSUBI = 0;              ///< SUBI miscibility flag
    int nParamCS = 0;               ///< Number of parameters
    int nCountSublatticeCS = 0;     ///< Sublattice count
    int nMaxSpeciesPhaseCS = 0;     ///< Max species per phase
    int nMagParamCS = 0;            ///< Magnetic parameters

    // 1D integer arrays
    Eigen::VectorXi nSpeciesPhaseCS;        ///< Species per phase
    Eigen::VectorXi nGibbsEqSpecies;        ///< Gibbs equations per species
    Eigen::VectorXi iPhaseCS;               ///< Phase indices
    Eigen::VectorXi iParticlesPerMoleCS;    ///< Particles per mole
    Eigen::VectorXi nParamPhaseCS;          ///< Parameters per phase
    Eigen::VectorXi iParamPassCS;           ///< Parameter pass flags
    Eigen::VectorXi nSublatticePhaseCS;     ///< Sublattices per phase
    Eigen::VectorXi iPhaseSublatticeCS;     ///< Sublattice ID per phase
    Eigen::VectorXi iMagParamPassCS;        ///< Magnetic param pass
    Eigen::VectorXi nMagParamPhaseCS;       ///< Magnetic params per phase
    Eigen::VectorXi iSUBIMixTypeCS;         ///< SUBI mix types
    Eigen::VectorXi nInterpolationOverrideCS; ///< Interpolation overrides

    // 2D integer matrices
    Eigen::MatrixXi iRegularParamCS;
    Eigen::MatrixXi nConstituentSublatticeCS;
    Eigen::MatrixXi nPairsSROCS;
    Eigen::MatrixXi iMagneticParamCS;
    Eigen::MatrixXi iSUBIParamDataCS;

    // 3D integer arrays
    std::vector<Eigen::MatrixXi> iInterpolationOverrideCS;
    std::vector<Eigen::MatrixXi> iConstituentSublatticeCS;
    std::vector<Eigen::MatrixXi> iPairIDCS;
    std::vector<Eigen::MatrixXi> iChemicalGroupCS;

    // 1D real arrays
    Eigen::VectorXd dAtomicMassCS;

    // 2D real matrices
    Eigen::MatrixXd dGibbsCoeffSpeciesTemp;  ///< Gibbs coefficients [species x coeffs]
    Eigen::MatrixXd dRegularParamCS;
    Eigen::MatrixXd dGibbsMagneticCS;
    Eigen::MatrixXd dMagneticParamCS;
    Eigen::MatrixXd dStoichSublatticeCS;
    Eigen::MatrixXd dStoichSpeciesCS;
    Eigen::MatrixXd dZetaSpeciesCS;
    Eigen::MatrixXd dStoichConstituentCS;
    Eigen::MatrixXd dQKTOParamsCS;

    // 3D real arrays
    std::vector<Eigen::MatrixXd> dSublatticeChargeCS;
    std::vector<Eigen::MatrixXd> dStoichPairsCS;
    std::vector<Eigen::MatrixXd> dConstituentCoefficientsCS;
    std::vector<Eigen::MatrixXd> dCoordinationNumberCS;

    // Character arrays
    std::vector<std::string> cElementNameCS;     ///< Element names (3 chars)
    std::vector<std::string> cSolnPhaseTypeCS;   ///< Phase types (8 chars)
    std::vector<std::string> cSolnPhaseNameCS;   ///< Phase names (25 chars)
    std::vector<std::string> cSpeciesNameCS;     ///< Species names (25 chars)
    std::vector<char> cRegularParamCS;           ///< Regular param flags

    // 2D/3D character arrays
    std::vector<std::vector<std::string>> cPairNameCS;
    std::vector<std::vector<std::vector<std::string>>> cConstituentNameSUBCS;

    // Constructor
    ParserState() = default;

    /// Allocate arrays for parsing
    void allocate(int nElements, int nSpecies, int nSolnPhases);

    /// Clear all parser state
    void clear();

    /// Check if phase type is supported
    static bool isPhaseTypeSupported(const std::string& type);

    /// Get phase type enum from string
    static Constants::PhaseType getPhaseType(const std::string& type);
};

} // namespace Thermochimica
