#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <array>
#include "../util/Constants.hpp"
#include "../util/Tolerances.hpp"

namespace Thermochimica {

/// Core thermodynamic state - converted from ModuleThermo.f90
/// Contains the main data structures for thermodynamic calculations
struct ThermoState {
    // System dimensions
    int nElements = 0;              ///< Number of elements in the system
    int nSpecies = 0;               ///< Number of species in the system
    int nParam = 0;                 ///< Number of mixing parameters
    int nMaxParam = 0;              ///< Maximum number of parameters allowed
    int nDummySpecies = 0;          ///< Number of dummy species from FactSage
    int nElemOrComp = 0;            ///< Elements or components
    int nMagParam = 0;              ///< Magnetic parameters

    // Phase counts
    int nConPhases = 0;             ///< Pure condensed phases at equilibrium
    int nSolnPhases = 0;            ///< Solution phases at equilibrium
    int nSolnPhasesSys = 0;         ///< Total solution phases in system
    int nConPhasesSys = 0;          ///< Total condensed phases in system
    int nChargedConstraints = 0;    ///< Electron constraints

    // Sublattice dimensions
    int nMaxSublatticeSys = 0;      ///< Maximum sublattices in system
    int nMaxConstituentSys = 0;     ///< Maximum constituents in system
    int nCountSublattice = 0;       ///< Number of sublattice phases

    // Scalar real values
    double dIdealConstant = Constants::kIdealGasConstant;  ///< Ideal gas constant
    double dNormalizeSum = Constants::kDefaultNormalize;    ///< Normalization sum
    double dNormalizeInput = Constants::kDefaultNormalize;  ///< Input normalization
    double dMassScale = 1.0;                                ///< Mass scaling factor
    double dTemperatureForLimits = 0.0;                     ///< Temperature for limits

    // Boolean flags
    bool lHeatCapacityCurrent = false;  ///< Heat capacity current flag

    // Tolerances
    Tolerances tolerances;

    // 1D integer arrays (allocatable in Fortran)
    Eigen::VectorXi iPhase;             ///< Phase type per species (0=pure, >0=soln, -1=dummy)
    Eigen::VectorXi nSpeciesPhase;      ///< Cumulative species count per phase (0-indexed)
    Eigen::VectorXi iParticlesPerMole;  ///< Particles per mole of constituent
    Eigen::VectorXi iSUBIMixType;       ///< SUBI mixing types
    Eigen::VectorXi iAssemblage;        ///< Phase indices in assemblage
    Eigen::VectorXi nParamPhase;        ///< Parameters per phase
    Eigen::VectorXi iElementSystem;     ///< Element system info
    Eigen::VectorXi iSpeciesPass;       ///< Species filter flags
    Eigen::VectorXi nMagParamPhase;     ///< Magnetic parameters per phase
    Eigen::VectorXi nSublatticePhase;   ///< Sublattices per phase
    Eigen::VectorXi iPhaseSublattice;   ///< Sublattice ID per phase
    Eigen::VectorXi iPhaseElectronID;   ///< Electron IDs
    Eigen::VectorXi nInterpolationOverride; ///< Interpolation overrides

    // 2D integer matrices
    Eigen::MatrixXi iRegularParam;          ///< Regular solution model info
    Eigen::MatrixXi iterHistoryLevel;       ///< Phase iteration history
    Eigen::MatrixXi nConstituentSublattice; ///< Constituents per sublattice
    Eigen::MatrixXi nPairsSRO;              ///< Short-range order pairs
    Eigen::MatrixXi iMagneticParam;         ///< Magnetic parameter data
    Eigen::MatrixXi iSUBLParamData;         ///< Sublattice parameter data

    // 3D integer arrays (stored as vector of matrices)
    std::vector<Eigen::MatrixXi> iInterpolationOverride;
    std::vector<Eigen::MatrixXi> iConstituentPass;
    std::vector<Eigen::MatrixXi> iConstituentSublattice;
    std::vector<Eigen::MatrixXi> iPairID;
    std::vector<Eigen::MatrixXi> iChemicalGroup;

    // 1D real arrays
    Eigen::VectorXd dStdGibbsEnergy;     ///< Standard Gibbs energy per species
    Eigen::VectorXd dGibbsSolnPhase;     ///< Gibbs energy per solution phase
    Eigen::VectorXd dMolesSpecies;       ///< Moles of each species
    Eigen::VectorXd dMagGibbsEnergy;     ///< Magnetic Gibbs energy
    Eigen::VectorXd dChemicalPotential;  ///< Chemical potential per species
    Eigen::VectorXd dExcessGibbsParam;   ///< Excess Gibbs parameters
    Eigen::VectorXd dLevel;              ///< Element potential adjustments
    Eigen::VectorXd dSpeciesTotalAtoms;  ///< Total atoms per species
    Eigen::VectorXd dElementPotential;   ///< Element potentials
    Eigen::VectorXd dMolesPhase;         ///< Moles of each phase
    Eigen::VectorXd dMolesElement;       ///< Moles of each element
    Eigen::VectorXd dMolFraction;        ///< Mole fractions
    Eigen::VectorXd dAtomicMass;         ///< Atomic masses

    // 2D real matrices
    Eigen::MatrixXd dAtomFractionSpecies;   ///< Atom fraction matrix
    Eigen::MatrixXd dStoichSublattice;      ///< Sublattice stoichiometry
    Eigen::MatrixXd dStoichSpecies;         ///< Stoichiometry matrix (nSpecies x nElements)
    Eigen::MatrixXd dStoichSpeciesUnFuzzed; ///< Unfuzzed stoichiometry
    Eigen::MatrixXd dQKTOParams;            ///< QKTO model parameters
    Eigen::MatrixXd dCoeffGibbsMagnetic;    ///< Magnetic Gibbs coefficients
    Eigen::MatrixXd dZetaSpecies;           ///< Zeta species parameters
    Eigen::MatrixXd dMagneticParam;         ///< Magnetic parameters

    // 3D real arrays (stored as vector of matrices)
    std::vector<Eigen::MatrixXd> dSiteFraction;          ///< Site fractions [phase][sublattice x constituent]
    std::vector<Eigen::MatrixXd> dCoordinationNumber;    ///< Coordination numbers
    std::vector<Eigen::MatrixXd> dSublatticeCharge;      ///< Sublattice charges
    std::vector<Eigen::MatrixXd> dStoichPairs;           ///< Stoichiometry pairs
    std::vector<Eigen::MatrixXd> dConstituentCoefficients; ///< Constituent coefficients

    // Character arrays (using std::string)
    std::vector<std::string> cElementName;      ///< Element names (max 12 chars)
    std::vector<std::string> cSpeciesName;      ///< Species names (max 30 chars)
    std::vector<std::string> cSolnPhaseType;    ///< Phase types (max 8 chars)
    std::vector<std::string> cSolnPhaseName;    ///< Phase names (max 25 chars)
    std::vector<char> cRegularParam;            ///< Regular param flags

    // 2D/3D character arrays
    std::vector<std::vector<std::string>> cPairName;           ///< Pair names [phase][pair]
    std::vector<std::vector<std::vector<std::string>>> cConstituentNameSUB; ///< Constituent names [phase][sublattice][constituent]

    // Constructor
    ThermoState() = default;

    /// Allocate arrays based on system dimensions
    void allocate(int numElements, int numSpecies, int numSolnPhases, int numParams);

    /// Deallocate all arrays
    void deallocate();

    /// Reset to initial state
    void reset();

    /// Get species index by name (returns -1 if not found)
    int getSpeciesIndex(const std::string& name) const;

    /// Get phase index by name (returns -1 if not found)
    int getPhaseIndex(const std::string& name) const;

    /// Get element index by name (returns -1 if not found)
    int getElementIndex(const std::string& name) const;
};

} // namespace Thermochimica
