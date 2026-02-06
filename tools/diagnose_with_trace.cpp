/// @file diagnose_with_trace.cpp
/// @brief Enhanced diagnostic with detailed solver trace
/// @details Adds debug output to understand why pure carbon produces 0 moles

#include "thermochimica/ThermoClass.hpp"
#include "thermochimica/Thermochimica.hpp"
#include <iostream>
#include <iomanip>

using namespace Thermochimica;

void testPureCarbonDetailed() {
    std::cout << "DETAILED TRACE: Pure Carbon System\n";
    std::cout << std::string(70, '=') << "\n\n";

    ThermoClass thermo;

    std::cout << "Step 1: Load database\n";
    int loadResult = thermo.loadDatabase("CO.dat");
    if (loadResult != 0) {
        std::cout << "❌ Failed to load database: " << thermo.getErrorMessage() << "\n";
        return;
    }
    std::cout << "✓ Database loaded\n\n";

    std::cout << "Step 2: Set inputs (C=1 mol, O=0 mol, T=900K, P=1 atm)\n";
    thermo.setStandardUnits();
    thermo.setTemperaturePressure(900.0, 1.0);
    thermo.setElementMass("C", 1.0);
    // No oxygen
    std::cout << "✓ Inputs set\n\n";

    // Get initial state
    const auto& ctx = thermo.getContext();
    const auto* state = ctx.thermo.get();

    std::cout << "Step 3: System info before calculation\n";
    std::cout << "  Total elements: " << state->nElements << "\n";
    std::cout << "  Total species: " << state->nSpecies << "\n";
    std::cout << "  Solution phases in database: " << state->nSolnPhasesSys << "\n";
    std::cout << "  Pure condensed species: "
              << (state->nSpecies - state->nSpeciesPhase[state->nSolnPhasesSys]) << "\n\n";

    std::cout << "Step 4: Calculate\n";
    int result = thermo.calculate();

    if (result != 0) {
        std::cout << "❌ Calculation failed: " << thermo.getErrorMessage() << "\n";
        std::cout << "   Error code: " << result << "\n";
    } else {
        std::cout << "✓ Calculation completed\n";
    }
    std::cout << "\n";

    std::cout << "Step 5: Results\n";
    std::cout << std::string(70, '-') << "\n";

    // System state
    std::cout << "System state:\n";
    std::cout << "  nConPhases (pure condensed): " << state->nConPhases << "\n";
    std::cout << "  nSolnPhases (solution): " << state->nSolnPhases << "\n";
    std::cout << "  Total phases: " << state->nConPhases + state->nSolnPhases << "\n\n";

    // Element mass balance
    std::cout << "Element mass balance:\n";
    std::cout << "  Input C: 1.0 mol\n";
    std::cout << "  Input O: 0.0 mol\n";

    double totalC = 0.0;
    double totalO = 0.0;

    // Find element indices
    int cIndex = -1, oIndex = -1;
    for (int j = 0; j < state->nElements; ++j) {
        std::string elemName = state->cElementName[j];
        if (elemName == "C") cIndex = j;
        if (elemName == "O") oIndex = j;
    }

    // Calculate total from phases
    for (int i = 0; i < state->nSpecies; ++i) {
        double moles = state->dMolesSpecies(i);
        if (moles > 1e-30) {
            if (cIndex >= 0) totalC += moles * state->dStoichSpecies(i, cIndex);
            if (oIndex >= 0) totalO += moles * state->dStoichSpecies(i, oIndex);
        }
    }

    std::cout << "  Output C: " << totalC << " mol\n";
    std::cout << "  Output O: " << totalO << " mol\n";
    std::cout << "  C balance error: " << std::abs(1.0 - totalC) << " mol\n\n";

    // Phase assemblage
    std::cout << "Phase assemblage:\n";

    // Solution phases
    for (int i = 0; i < state->nSolnPhases; ++i) {
        int idx = state->nElements - state->nSolnPhases + i;
        int phaseID = state->iAssemblage[idx];
        if (phaseID < 0) {
            int solnIdx = -phaseID - 1;
            double moles = state->dMolesPhase[idx];
            std::cout << "  [" << i << "] Solution: " << state->cSolnPhaseName[solnIdx]
                      << " = " << moles << " mol\n";

            // Show species in this phase
            int iFirst = (solnIdx > 0) ? state->nSpeciesPhase[solnIdx] : 0;
            int iLast = state->nSpeciesPhase[solnIdx + 1];
            for (int j = iFirst; j < iLast; ++j) {
                if (state->dMolesSpecies[j] > 1e-30) {
                    std::cout << "      - " << state->cSpeciesName[j]
                              << ": " << state->dMolesSpecies[j] << " mol\n";
                }
            }
        }
    }

    // Pure condensed phases
    for (int i = 0; i < state->nConPhases; ++i) {
        int speciesIdx = state->iAssemblage[i] - 1;
        if (speciesIdx >= 0 && speciesIdx < state->nSpecies) {
            double moles = state->dMolesPhase[i];
            std::cout << "  [" << i << "] Pure: " << state->cSpeciesName[speciesIdx]
                      << " = " << moles << " mol\n";
        }
    }

    if (state->nConPhases == 0 && state->nSolnPhases == 0) {
        std::cout << "  ⚠️ NO PHASES IN ASSEMBLAGE!\n";
    }
    std::cout << "\n";

    // Chemical potentials
    auto [muC, infoC] = thermo.getOutputChemPot("C");
    std::cout << "Chemical potentials:\n";
    std::cout << "  μ_C = " << muC << " J/mol\n\n";

    // Gibbs energy
    std::cout << "Thermodynamics:\n";
    std::cout << "  G_total = " << thermo.getGibbsEnergy() << " J\n\n";

    // Check driving forces manually
    std::cout << "Driving force check (manual):\n";

    // Find graphite species
    for (int i = state->nSpeciesPhase[state->nSolnPhasesSys]; i < state->nSpecies; ++i) {
        std::string name = state->cSpeciesName[i];
        if (name.find("C") != std::string::npos || name.find("graphite") != std::string::npos) {
            double gStd = state->dStdGibbsEnergy[i];
            double muCalc = muC;  // For pure carbon, μ_calc = μ_C
            double drivingForce = gStd - muCalc;

            std::cout << "  Species: " << name << "\n";
            std::cout << "    G_std = " << gStd << " J/mol\n";
            std::cout << "    μ_calc = " << muCalc << " J/mol\n";
            std::cout << "    DF = G_std - μ_calc = " << drivingForce << " J/mol\n";
            std::cout << "    Stable? " << (drivingForce < 0 ? "YES" : "NO") << "\n\n";
        }
    }

    std::cout << "ANALYSIS:\n";
    std::cout << std::string(70, '-') << "\n";

    if (std::abs(1.0 - totalC) > 0.01) {
        std::cout << "❌ MASS CONSERVATION VIOLATED!\n";
        std::cout << "   Input: 1.0 mol C\n";
        std::cout << "   Output: " << totalC << " mol C\n";
        std::cout << "   This is thermodynamically impossible!\n\n";
    }

    if (state->nConPhases == 0 && state->nSolnPhases == 0) {
        std::cout << "❌ NO PHASES STABLE!\n";
        std::cout << "   This violates the phase rule.\n";
        std::cout << "   For a 1-component system at fixed T,P, at least 1 phase must exist.\n\n";
    }

    if (muC > 0) {
        std::cout << "⚠️ POSITIVE CHEMICAL POTENTIAL!\n";
        std::cout << "   μ_C > 0 indicates graphite should be stable.\n";
        std::cout << "   But it's not in the assemblage.\n\n";
    }

    std::cout << "\n" << std::string(70, '=') << "\n";
}

int main() {
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║         ENHANCED DIAGNOSTIC WITH DETAILED TRACE                  ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";

    testPureCarbonDetailed();

    return 0;
}
