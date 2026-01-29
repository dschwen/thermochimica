#include "thermochimica/Thermochimica.hpp"
#include <iostream>
#include <iomanip>

namespace Thermochimica {

void printResultsDetailed(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;
    auto& io = *ctx.io;

    std::cout << "\n";
    std::cout << "========================================\n";
    std::cout << "        THERMOCHIMICA RESULTS\n";
    std::cout << "========================================\n";

    std::cout << "\nConditions:\n";
    std::cout << "  Temperature: " << std::fixed << std::setprecision(2)
              << io.dTemperature << " K\n";
    std::cout << "  Pressure:    " << io.dPressure << " atm\n";

    std::cout << "\nSystem Gibbs Energy: " << std::scientific << std::setprecision(6)
              << io.dGibbsEnergySys << " J\n";

    if (thermo.nSolnPhases > 0) {
        std::cout << "\nStable Solution Phases:\n";
        for (int i = 0; i < thermo.nSolnPhases; ++i) {
            int idx = thermo.nElements - thermo.nSolnPhases + i;
            int phaseIdx = -thermo.iAssemblage(idx) - 1;
            if (phaseIdx >= 0) {
                std::cout << "  " << std::setw(20) << std::left
                          << thermo.cSolnPhaseName[phaseIdx]
                          << ": " << std::scientific << std::setprecision(4)
                          << thermo.dMolesPhase(idx) << " mol\n";
            }
        }
    }

    if (thermo.nConPhases > 0) {
        std::cout << "\nStable Pure Condensed Phases:\n";
        for (int i = 0; i < thermo.nConPhases; ++i) {
            int speciesIdx = thermo.iAssemblage(i) - 1;
            if (speciesIdx >= 0) {
                std::cout << "  " << std::setw(20) << std::left
                          << thermo.cSpeciesName[speciesIdx]
                          << ": " << std::scientific << std::setprecision(4)
                          << thermo.dMolesPhase(i) << " mol\n";
            }
        }
    }

    std::cout << "\n========================================\n\n";
}

} // namespace Thermochimica
