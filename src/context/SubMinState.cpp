#include "thermochimica/context/SubMinState.hpp"

namespace Thermochimica {

void SubMinState::allocate(int nSpecies) {
    nVar = nSpecies;

    iHessian.resize(nSpecies);
    iHessian.setZero();

    dChemicalPotentialStar.resize(nSpecies);
    dChemicalPotentialStar.setZero();

    dRHS.resize(nSpecies);
    dRHS.setZero();

    dHessian.resize(nSpecies, nSpecies);
    dHessian.setZero();
}

void SubMinState::reset() {
    nVar = 0;
    iFirst = 0;
    iLast = 0;
    iSolnPhaseIndexOther = 0;

    dDrivingForce = 0.0;
    dDrivingForceLast = 0.0;
    dSubMinFunctionNorm = 0.0;

    lSubMinConverged = false;

    if (iHessian.size() > 0) {
        iHessian.setZero();
    }
    if (dChemicalPotentialStar.size() > 0) {
        dChemicalPotentialStar.setZero();
    }
    if (dRHS.size() > 0) {
        dRHS.setZero();
    }
    if (dHessian.size() > 0) {
        dHessian.setZero();
    }
}

void SubMinState::initForPhase(int first, int last) {
    iFirst = first;
    iLast = last;
    nVar = last - first + 1;

    dDrivingForce = 0.0;
    dDrivingForceLast = 0.0;
    dSubMinFunctionNorm = 0.0;
    lSubMinConverged = false;

    // Resize working arrays if needed
    if (dRHS.size() < nVar) {
        dRHS.resize(nVar);
    }
    if (dHessian.rows() < nVar || dHessian.cols() < nVar) {
        dHessian.resize(nVar, nVar);
    }
    if (dChemicalPotentialStar.size() < nVar) {
        dChemicalPotentialStar.resize(nVar);
    }

    dRHS.head(nVar).setZero();
    dHessian.topLeftCorner(nVar, nVar).setZero();
    dChemicalPotentialStar.head(nVar).setZero();
}

} // namespace Thermochimica
