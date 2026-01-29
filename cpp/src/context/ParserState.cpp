#include "thermochimica/context/ParserState.hpp"
#include <algorithm>

namespace Thermochimica {

void ParserState::allocate(int nElements, int nSpecies, int nSolnPhases) {
    nElementsCS = nElements;
    nSpeciesCS = nSpecies;
    nSolnPhasesSysCS = nSolnPhases;

    // 1D integer arrays
    nSpeciesPhaseCS.resize(nSolnPhases + 1);
    nSpeciesPhaseCS.setZero();

    nGibbsEqSpecies.resize(nSpecies);
    nGibbsEqSpecies.setOnes();  // Default 1 equation per species

    iPhaseCS.resize(nSpecies);
    iPhaseCS.setZero();

    iParticlesPerMoleCS.resize(nSpecies);
    iParticlesPerMoleCS.setOnes();

    nParamPhaseCS.resize(nSolnPhases + 1);
    nParamPhaseCS.setZero();

    iParamPassCS.resize(nSolnPhases);
    iParamPassCS.setZero();

    nSublatticePhaseCS.resize(nSolnPhases);
    nSublatticePhaseCS.setZero();

    iPhaseSublatticeCS.resize(nSolnPhases);
    iPhaseSublatticeCS.setZero();

    iMagParamPassCS.resize(nSolnPhases);
    iMagParamPassCS.setZero();

    nMagParamPhaseCS.resize(nSolnPhases + 1);
    nMagParamPhaseCS.setZero();

    iSUBIMixTypeCS.resize(nSolnPhases);
    iSUBIMixTypeCS.setZero();

    nInterpolationOverrideCS.resize(nSolnPhases + 1);  // +1 for 1-based indexing
    nInterpolationOverrideCS.setZero();

    // 1D real arrays
    dAtomicMassCS.resize(nElements);
    dAtomicMassCS.setZero();

    // 2D real matrices
    dGibbsCoeffSpeciesTemp.resize(nSpecies, kNumGibbsCoeff * kMaxGibbsEqs);
    dGibbsCoeffSpeciesTemp.setZero();

    dStoichSpeciesCS.resize(nSpecies, nElements);
    dStoichSpeciesCS.setZero();

    // Character arrays
    cElementNameCS.resize(nElements);
    cSolnPhaseTypeCS.resize(nSolnPhases);
    cSolnPhaseNameCS.resize(nSolnPhases);
    cSpeciesNameCS.resize(nSpecies);
}

void ParserState::clear() {
    nElementsCS = 0;
    nSpeciesCS = 0;
    nSolnPhasesSysCS = 0;
    INFO = 0;
    iMiscSUBI = 0;
    nParamCS = 0;
    nCountSublatticeCS = 0;
    nMaxSpeciesPhaseCS = 0;
    nMagParamCS = 0;

    // Clear all Eigen vectors/matrices (automatic)
    nSpeciesPhaseCS.resize(0);
    nGibbsEqSpecies.resize(0);
    iPhaseCS.resize(0);
    iParticlesPerMoleCS.resize(0);
    nParamPhaseCS.resize(0);
    iParamPassCS.resize(0);
    nSublatticePhaseCS.resize(0);
    iPhaseSublatticeCS.resize(0);

    dAtomicMassCS.resize(0);
    dGibbsCoeffSpeciesTemp.resize(0, 0);
    dStoichSpeciesCS.resize(0, 0);

    // Clear string vectors
    cElementNameCS.clear();
    cSolnPhaseTypeCS.clear();
    cSolnPhaseNameCS.clear();
    cSpeciesNameCS.clear();
    cRegularParamCS.clear();
    cPairNameCS.clear();
    cConstituentNameSUBCS.clear();

    // Clear 3D arrays
    iInterpolationOverrideCS.clear();
    iConstituentSublatticeCS.clear();
    iPairIDCS.clear();
    iChemicalGroupCS.clear();
    dSublatticeChargeCS.clear();
    dStoichPairsCS.clear();
    dConstituentCoefficientsCS.clear();
    dCoordinationNumberCS.clear();
}

bool ParserState::isPhaseTypeSupported(const std::string& type) {
    std::string upperType = type;
    std::transform(upperType.begin(), upperType.end(), upperType.begin(), ::toupper);

    // Trim to 8 characters and compare
    if (upperType.size() > 8) {
        upperType = upperType.substr(0, 8);
    }

    for (int i = 0; i < kNumSolnTypeSupport; ++i) {
        std::string supported = kSolnPhaseTypeSupport[i];
        if (upperType == supported || upperType.find(supported) == 0) {
            return true;
        }
    }
    return false;
}

Constants::PhaseType ParserState::getPhaseType(const std::string& type) {
    std::string upperType = type;
    std::transform(upperType.begin(), upperType.end(), upperType.begin(), ::toupper);

    // Trim whitespace
    while (!upperType.empty() && upperType.back() == ' ') {
        upperType.pop_back();
    }

    if (upperType == "IDMX") return Constants::PhaseType::IDMX;
    if (upperType == "QKTO") return Constants::PhaseType::QKTO;
    if (upperType == "SUBL") return Constants::PhaseType::SUBL;
    if (upperType == "RKMP") return Constants::PhaseType::RKMP;
    if (upperType == "RKMPM") return Constants::PhaseType::RKMPM;
    if (upperType == "SUBLM") return Constants::PhaseType::SUBLM;
    if (upperType == "SUBG") return Constants::PhaseType::SUBG;
    if (upperType == "SUBQ") return Constants::PhaseType::SUBQ;
    if (upperType == "SUBI") return Constants::PhaseType::SUBI;
    if (upperType == "SUBM") return Constants::PhaseType::SUBM;

    return Constants::PhaseType::Unknown;
}

} // namespace Thermochimica
