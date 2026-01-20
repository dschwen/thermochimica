#include "thermochimica/context/GEMState.hpp"

namespace Thermochimica {

void GEMState::allocate(int nElements, int nSolnPhasesSys, int nSpecies) {
    // Iteration history
    iterHistory.resize(nElements, Constants::kIterGlobalMax);
    iterHistory.setZero();

    // Working vectors
    dSumMolFractionSoln.resize(nSolnPhasesSys);
    dSumMolFractionSoln.setZero();

    dMolesPhaseLast.resize(nElements);
    dMolesPhaseLast.setZero();

    // Update vector size = nElements + nSolnPhases + nConPhases (max = 2*nElements)
    dUpdateVar.resize(2 * nElements);
    dUpdateVar.setZero();

    dDrivingForceSoln.resize(nSolnPhasesSys);
    dDrivingForceSoln.setZero();

    dPartialExcessGibbs.resize(nSpecies);
    dPartialExcessGibbs.setZero();

    dPartialExcessGibbsLast.resize(nSpecies);
    dPartialExcessGibbsLast.setZero();

    // Effective stoichiometry matrix
    dEffStoichSolnPhase.resize(nSolnPhasesSys, nElements);
    dEffStoichSolnPhase.setZero();

    // Phase status arrays
    lSolnPhases.resize(nSolnPhasesSys, false);
    lMiscibility.resize(nSolnPhasesSys, false);
}

void GEMState::reset() {
    iterLast = 0;
    iterStep = 0;
    iterRevert = 0;
    iterGlobal = 0;
    iterLastCon = 0;
    iterLastSoln = 0;
    iterSwap = 0;
    iterLastMiscGapCheck = 0;

    iConPhaseLast = 0;
    iSolnPhaseLast = 0;
    iSolnSwap = 0;
    iPureConSwap = 0;

    dGEMFunctionNorm = 0.0;
    dGEMFunctionNormLast = 0.0;
    dMaxSpeciesChange = 0.0;
    dMinGibbs = 0.0;

    lDebugMode = false;
    lRevertSystem = false;
    lConverged = false;

    // Reset arrays if allocated
    if (iterHistory.size() > 0) {
        iterHistory.setZero();
    }
    if (dSumMolFractionSoln.size() > 0) {
        dSumMolFractionSoln.setZero();
    }
    if (dMolesPhaseLast.size() > 0) {
        dMolesPhaseLast.setZero();
    }
    if (dUpdateVar.size() > 0) {
        dUpdateVar.setZero();
    }
    if (dDrivingForceSoln.size() > 0) {
        dDrivingForceSoln.setZero();
    }
    if (dPartialExcessGibbs.size() > 0) {
        dPartialExcessGibbs.setZero();
    }
    if (dPartialExcessGibbsLast.size() > 0) {
        dPartialExcessGibbsLast.setZero();
    }
    if (dEffStoichSolnPhase.size() > 0) {
        dEffStoichSolnPhase.setZero();
    }

    std::fill(lSolnPhases.begin(), lSolnPhases.end(), false);
    std::fill(lMiscibility.begin(), lMiscibility.end(), false);
}

void GEMState::initIteration() {
    // Called at start of each GEM iteration
    dGEMFunctionNormLast = dGEMFunctionNorm;

    if (dMolesPhaseLast.size() > 0) {
        // Save current moles for comparison
        // This will be set from thermo->dMolesPhase in the actual solver
    }

    if (dPartialExcessGibbsLast.size() > 0 && dPartialExcessGibbs.size() > 0) {
        dPartialExcessGibbsLast = dPartialExcessGibbs;
    }
}

} // namespace Thermochimica
