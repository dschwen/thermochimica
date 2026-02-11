#include "thermochimica/Thermochimica.hpp"
#include <cmath>

namespace Thermochimica {

// Forward declaration of internal function
void computeHeatCapacity(ThermoContext& ctx);

void postProcess(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;
    auto& io = *ctx.io;

    // Calculate system Gibbs energy
    double G = 0.0;

    double RT = Constants::kIdealGasConstant * io.dTemperature;

    // Contribution from solution phases
    for (int iPhase = 0; iPhase < thermo.nSolnPhases; ++iPhase) {
        int phaseIdx = -thermo.iAssemblage(thermo.nElements - thermo.nSolnPhases + iPhase) - 1;
        if (phaseIdx < 0) continue;

        int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
        int iLast = thermo.nSpeciesPhase(phaseIdx + 1);

        for (int i = iFirst; i < iLast; ++i) {
            double n = thermo.dMolesSpecies(i);
            double mu = thermo.dChemicalPotential(i);
            if (n > 1e-100) {
                G += n * mu;
            }
        }
    }

    // Contribution from pure condensed phases
    for (int i = 0; i < thermo.nConPhases; ++i) {
        int speciesIdx = thermo.iAssemblage(i) - 1;
        if (speciesIdx >= 0 && speciesIdx < thermo.nSpecies) {
            G += thermo.dMolesPhase(i) * thermo.dChemicalPotential(speciesIdx);
        }
    }

    io.dGibbsEnergySys = G * RT;

    // Prepare output arrays
    io.nSolnPhasesOut = thermo.nSolnPhases;
    io.nPureConPhaseOut = thermo.nConPhases;

    io.cSolnPhaseNameOut.clear();
    io.dSolnPhaseMolesOut.clear();

    for (int i = 0; i < thermo.nSolnPhases; ++i) {
        int idx = thermo.nElements - thermo.nSolnPhases + i;
        int phaseIdx = -thermo.iAssemblage(idx) - 1;
        if (phaseIdx >= 0 && phaseIdx < static_cast<int>(thermo.cSolnPhaseName.size())) {
            io.cSolnPhaseNameOut.push_back(thermo.cSolnPhaseName[phaseIdx]);
            io.dSolnPhaseMolesOut.push_back(thermo.dMolesPhase(idx));
        }
    }

    io.cPureConPhaseNameOut.clear();
    io.dPureConPhaseMolesOut.clear();

    for (int i = 0; i < thermo.nConPhases; ++i) {
        int speciesIdx = thermo.iAssemblage(i) - 1;
        if (speciesIdx >= 0 && speciesIdx < static_cast<int>(thermo.cSpeciesName.size())) {
            io.cPureConPhaseNameOut.push_back(thermo.cSpeciesName[speciesIdx]);
            io.dPureConPhaseMolesOut.push_back(thermo.dMolesPhase(i));
        }
    }

    // Compute heat capacity, entropy, and enthalpy if requested
    computeHeatCapacity(ctx);
}

} // namespace Thermochimica
