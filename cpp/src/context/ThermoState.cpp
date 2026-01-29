#include "thermochimica/context/ThermoState.hpp"
#include <algorithm>

namespace Thermochimica {

void ThermoState::allocate(int numElements, int numSpecies, int numSolnPhases, int numParams) {
    nElements = numElements;
    nSpecies = numSpecies;
    nSolnPhasesSys = numSolnPhases;
    nParam = numParams;

    // 1D integer arrays
    iPhase.resize(nSpecies);
    iPhase.setZero();

    nSpeciesPhase.resize(nSolnPhasesSys + 1);
    nSpeciesPhase.setZero();

    iParticlesPerMole.resize(nSpecies);
    iParticlesPerMole.setOnes();  // Default is 1

    iAssemblage.resize(nElements);
    iAssemblage.setZero();

    nParamPhase.resize(nSolnPhasesSys + 1);
    nParamPhase.setZero();

    iElementSystem.resize(nElements);
    iElementSystem.setZero();

    iSpeciesPass.resize(nSpecies);
    iSpeciesPass.setOnes();  // All species pass by default

    nMagParamPhase.resize(nSolnPhasesSys + 1);
    nMagParamPhase.setZero();

    nSublatticePhase.resize(nSolnPhasesSys);
    nSublatticePhase.setZero();

    iPhaseSublattice.resize(nSolnPhasesSys);
    iPhaseSublattice.setZero();

    iPhaseElectronID.resize(nSolnPhasesSys);
    iPhaseElectronID.setZero();

    iSUBIMixType.resize(nSolnPhasesSys);
    iSUBIMixType.setZero();

    nInterpolationOverride.resize(nSolnPhasesSys + 1);  // +1 for 1-based indexing
    nInterpolationOverride.setZero();

    // 2D integer matrices
    if (numParams > 0) {
        iRegularParam.resize(numParams, Constants::kMaxParams + 2);
        iRegularParam.setZero();
    }

    iterHistoryLevel.resize(nElements, 1000);  // History size
    iterHistoryLevel.setZero();

    // 1D real arrays
    dStdGibbsEnergy.resize(nSpecies);
    dStdGibbsEnergy.setZero();

    dGibbsSolnPhase.resize(nSolnPhasesSys);
    dGibbsSolnPhase.setZero();

    dMolesSpecies.resize(nSpecies);
    dMolesSpecies.setZero();

    dMagGibbsEnergy.resize(nSpecies);
    dMagGibbsEnergy.setZero();

    dChemicalPotential.resize(nSpecies);
    dChemicalPotential.setZero();

    if (numParams > 0) {
        dExcessGibbsParam.resize(numParams);
        dExcessGibbsParam.setZero();
    }

    dLevel.resize(nElements);
    dLevel.setZero();

    dSpeciesTotalAtoms.resize(nSpecies);
    dSpeciesTotalAtoms.setZero();

    dElementPotential.resize(nElements);
    dElementPotential.setZero();

    dMolesPhase.resize(nElements);
    dMolesPhase.setZero();

    dMolesElement.resize(nElements);
    dMolesElement.setZero();

    dMolFraction.resize(nSpecies);
    dMolFraction.setZero();

    dAtomicMass.resize(nElements);
    dAtomicMass.setZero();

    // 2D real matrices
    dStoichSpecies.resize(nSpecies, nElements);
    dStoichSpecies.setZero();

    dStoichSpeciesUnFuzzed.resize(nSpecies, nElements);
    dStoichSpeciesUnFuzzed.setZero();

    // QKTO parameters (3 columns: K, p, group)
    dQKTOParams.resize(nSpecies, 3);
    dQKTOParams.setZero();

    dAtomFractionSpecies.resize(nSpecies, nElements);
    dAtomFractionSpecies.setZero();

    // Character arrays
    cElementName.resize(nElements);
    cSpeciesName.resize(nSpecies);
    cSolnPhaseType.resize(nSolnPhasesSys);
    cSolnPhaseName.resize(nSolnPhasesSys);

    if (numParams > 0) {
        cRegularParam.resize(numParams, ' ');
    }

    // Initialize tolerances
    tolerances.initDefaults();
}

void ThermoState::deallocate() {
    // Eigen vectors/matrices automatically deallocate
    // Just reset dimensions
    nElements = 0;
    nSpecies = 0;
    nSolnPhasesSys = 0;
    nParam = 0;

    // Clear string vectors
    cElementName.clear();
    cSpeciesName.clear();
    cSolnPhaseType.clear();
    cSolnPhaseName.clear();
    cRegularParam.clear();
    cPairName.clear();
    cConstituentNameSUB.clear();
}

void ThermoState::reset() {
    deallocate();

    // Reset scalars
    nConPhases = 0;
    nSolnPhases = 0;
    nConPhasesSys = 0;
    nChargedConstraints = 0;
    nDummySpecies = 0;
    nMaxSublatticeSys = 0;
    nMaxConstituentSys = 0;
    nCountSublattice = 0;
    nMagParam = 0;
    nElemOrComp = 0;
    nMaxParam = 0;

    dIdealConstant = Constants::kIdealGasConstant;
    dNormalizeSum = Constants::kDefaultNormalize;
    dNormalizeInput = Constants::kDefaultNormalize;
    dMassScale = 1.0;
    dTemperatureForLimits = 0.0;
    lHeatCapacityCurrent = false;

    tolerances.initDefaults();
}

int ThermoState::getSpeciesIndex(const std::string& name) const {
    for (int i = 0; i < static_cast<int>(cSpeciesName.size()); ++i) {
        if (cSpeciesName[i] == name) {
            return i;
        }
    }
    return -1;
}

int ThermoState::getPhaseIndex(const std::string& name) const {
    for (int i = 0; i < static_cast<int>(cSolnPhaseName.size()); ++i) {
        if (cSolnPhaseName[i] == name) {
            return i;
        }
    }
    return -1;
}

int ThermoState::getElementIndex(const std::string& name) const {
    for (int i = 0; i < static_cast<int>(cElementName.size()); ++i) {
        if (cElementName[i] == name) {
            return i;
        }
    }
    return -1;
}

} // namespace Thermochimica
