#include "thermochimica/context/CTZState.hpp"

namespace Thermochimica {

void CTZState::allocate(int maxAssemblages, int maxElements, int maxPhases, int maxSpecies) {
    nMaxAssemblages = maxAssemblages;
    nMaxElements = maxElements;
    nAssemblages = 0;

    assemblageHistory.resize(maxAssemblages, maxElements);
    assemblageHistory.setZero();

    elementHistory.resize(maxAssemblages);
    for (auto& row : elementHistory) {
        row.resize(maxElements);
    }

    assemblageTlimits.resize(maxAssemblages, 2);
    assemblageTlimits.setZero();

    // 4D stoichiometry history: [assemblage][phase][species x element]
    stoichHistory.resize(maxAssemblages);
    for (auto& phaseVec : stoichHistory) {
        phaseVec.resize(maxPhases);
        for (auto& mat : phaseVec) {
            mat.resize(maxSpecies, maxElements);
            mat.setZero();
        }
    }

    lCtzInit = true;
}

void CTZState::init(int maxAssemblages, int maxElements) {
    // Simplified init with default phase/species dimensions
    allocate(maxAssemblages, maxElements, maxElements, 100);
}

void CTZState::reset() {
    nAssemblages = 0;
    maxNorm = 0.0;
    tRange = 0.0;

    if (assemblageHistory.size() > 0) {
        assemblageHistory.setZero();
    }
    if (assemblageTlimits.size() > 0) {
        assemblageTlimits.setZero();
    }

    for (auto& row : elementHistory) {
        std::fill(row.begin(), row.end(), std::string());
    }

    for (auto& phaseVec : stoichHistory) {
        for (auto& mat : phaseVec) {
            mat.setZero();
        }
    }

    // Keep lCtzInit true if already initialized
}

void CTZState::addAssemblage(const Eigen::VectorXi& assemblage,
                             const std::vector<std::string>& elements,
                             double tMin, double tMax) {
    if (nAssemblages >= nMaxAssemblages) {
        // Could implement circular buffer or resize
        return;
    }

    int idx = nAssemblages;

    // Store assemblage
    int nElem = std::min(static_cast<int>(assemblage.size()), nMaxElements);
    for (int i = 0; i < nElem; ++i) {
        assemblageHistory(idx, i) = assemblage(i);
    }

    // Store element names
    int nNames = std::min(static_cast<int>(elements.size()), nMaxElements);
    for (int i = 0; i < nNames; ++i) {
        elementHistory[idx][i] = elements[i];
    }

    // Store temperature limits
    assemblageTlimits(idx, 0) = tMin;
    assemblageTlimits(idx, 1) = tMax;

    ++nAssemblages;
}

} // namespace Thermochimica
