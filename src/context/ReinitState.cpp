#include "thermochimica/context/ReinitState.hpp"

namespace Thermochimica {

void ReinitState::allocate(int nElements, int nSpecies) {
    iAssemblage_Old.resize(nElements);
    iAssemblage_Old.setZero();

    iElementsUsed_Old.fill(0);

    dChemicalPotential_Old.resize(nSpecies);
    dChemicalPotential_Old.setZero();

    dMolesPhase_Old.resize(nElements);
    dMolesPhase_Old.setZero();

    dElementPotential_Old.resize(nElements);
    dElementPotential_Old.setZero();

    dMolFraction_Old.resize(nSpecies);
    dMolFraction_Old.setZero();
}

void ReinitState::saveState(const Eigen::VectorXi& assemblage,
                            const Eigen::VectorXd& chemPot,
                            const Eigen::VectorXd& molesPhase,
                            const Eigen::VectorXd& elemPot,
                            const Eigen::VectorXd& molFrac) {
    // Resize if necessary
    if (iAssemblage_Old.size() != assemblage.size()) {
        iAssemblage_Old.resize(assemblage.size());
    }
    if (dChemicalPotential_Old.size() != chemPot.size()) {
        dChemicalPotential_Old.resize(chemPot.size());
    }
    if (dMolesPhase_Old.size() != molesPhase.size()) {
        dMolesPhase_Old.resize(molesPhase.size());
    }
    if (dElementPotential_Old.size() != elemPot.size()) {
        dElementPotential_Old.resize(elemPot.size());
    }
    if (dMolFraction_Old.size() != molFrac.size()) {
        dMolFraction_Old.resize(molFrac.size());
    }

    iAssemblage_Old = assemblage;
    dChemicalPotential_Old = chemPot;
    dMolesPhase_Old = molesPhase;
    dElementPotential_Old = elemPot;
    dMolFraction_Old = molFrac;
}

bool ReinitState::isValid() const {
    // Check if we have valid saved state
    if (iAssemblage_Old.size() == 0) {
        return false;
    }

    // Check if any non-zero values exist
    bool hasData = false;
    for (int i = 0; i < iAssemblage_Old.size(); ++i) {
        if (iAssemblage_Old(i) != 0) {
            hasData = true;
            break;
        }
    }

    return hasData;
}

void ReinitState::clear() {
    if (iAssemblage_Old.size() > 0) {
        iAssemblage_Old.setZero();
    }

    iElementsUsed_Old.fill(0);

    if (dChemicalPotential_Old.size() > 0) {
        dChemicalPotential_Old.setZero();
    }
    if (dMolesPhase_Old.size() > 0) {
        dMolesPhase_Old.setZero();
    }
    if (dElementPotential_Old.size() > 0) {
        dElementPotential_Old.setZero();
    }
    if (dMolFraction_Old.size() > 0) {
        dMolFraction_Old.setZero();
    }
}

} // namespace Thermochimica
