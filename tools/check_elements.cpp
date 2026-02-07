/// @file check_elements.cpp
/// @brief Check element ordering in database

#include "thermochimica/ThermoClass.hpp"
#include <iostream>

using namespace Thermochimica;

int main() {
    ThermoClass thermo;
    thermo.loadDatabase("CO.dat");

    const auto& ctx = thermo.getContext();
    const auto* state = ctx.thermo.get();

    std::cout << "Elements in database:\n";
    for (int i = 0; i < state->nElements; ++i) {
        std::cout << "  [" << i << "] " << state->cElementName[i] << "\n";
    }

    std::cout << "\nAll pure condensed species:\n";
    for (int i = state->nSpeciesPhase[state->nSolnPhasesSys]; i < state->nSpecies; ++i) {
        std::cout << "  [" << i << "] " << state->cSpeciesName[i];
        std::cout << " (";
        for (int j = 0; j < state->nElements; ++j) {
            double stoich = state->dStoichSpecies(i, j);
            if (stoich > 0) {
                std::cout << state->cElementName[j] << ":" << stoich << " ";
            }
        }
        std::cout << ")\n";
    }

    return 0;
}
