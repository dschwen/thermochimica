/// @file IdealMixingModel.cpp
/// @brief Implementation of ideal mixing thermodynamic model

#include "thermochimica/models/IdealMixingModel.hpp"
#include "thermochimica/context/ThermoState.hpp"
#include "thermochimica/context/GEMState.hpp"
#include <cmath>

namespace Thermochimica {

void IdealMixingModel::computeExcessGibbs(ThermoState& thermo,
                                          GEMState& gem,
                                          int phaseIndex) {
    // For ideal mixing, excess Gibbs energy is zero
    thermo.dGibbsSolnPhase(phaseIndex) = 0.0;

    // Compute ideal mixing contribution to chemical potential
    // mu_i = mu_i^0 + RT*ln(x_i)

    int iFirst = (phaseIndex > 0) ? thermo.nSpeciesPhase(phaseIndex) : 0;
    int iLast = thermo.nSpeciesPhase(phaseIndex + 1);

    for (int i = iFirst; i < iLast; ++i) {
        double x = thermo.dMolFraction(i);
        if (x > 1e-100) {
            // Add RT*ln(x) contribution (already normalized by RT)
            thermo.dChemicalPotential(i) += std::log(x);
        }
    }
}

bool IdealMixingModel::canHandle(const ThermoState& thermo, int phaseIndex) const {
    if (phaseIndex < 0 || phaseIndex >= static_cast<int>(thermo.iSolnPhaseType.size())) {
        return false;
    }
    return thermo.iSolnPhaseType[phaseIndex] == Constants::PhaseType::IDMX;
}

} // namespace Thermochimica
