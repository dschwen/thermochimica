/// @file diagnose_graphite.cpp
/// @brief Diagnostic tests for graphite phase stability issue
/// @details Runs Tests A-D from PHASE1_DIAGNOSTICS.md

#include "thermochimica/ThermoClass.hpp"
#include "thermochimica/Thermochimica.hpp"
#include <iostream>
#include <iomanip>

using namespace Thermochimica;

void printSeparator() {
    std::cout << "\n" << std::string(70, '=') << "\n\n";
}

void testA_CheckGraphiteInDatabase() {
    std::cout << "TEST A: Check if Graphite is in Database\n";
    std::cout << std::string(70, '-') << "\n";

    ThermoClass thermo;
    int loadResult = thermo.loadDatabase("CO.dat");

    if (loadResult != 0) {
        std::cout << "❌ Failed to load database: " << thermo.getErrorMessage() << "\n";
        return;
    }

    std::cout << "✓ Database loaded successfully\n\n";

    // Try different phase name variations
    std::vector<std::string> phaseNames = {
        "C_Graphite(s)",
        "C_Graphite",
        "Graphite",
        "C(gr)",
        "C_gr"
    };

    // First do a minimal calculation to initialize
    thermo.setStandardUnits();
    thermo.setTemperaturePressure(900.0, 1.0);
    thermo.setElementMass("C", 1.0);
    thermo.setElementMass("O", 1.0);
    thermo.calculate();

    std::cout << "Testing phase name variations:\n";
    for (const auto& name : phaseNames) {
        auto [moles, info] = thermo.getMolesPhase(name);
        std::cout << "  " << std::setw(20) << std::left << name << ": ";
        if (info == 0) {
            std::cout << "Found! Moles = " << moles << "\n";
        } else {
            std::cout << "Not found (error " << info << ")\n";
        }
    }

    std::cout << "\n✓ Test A complete\n";
}

void testB_CheckChemicalPotential() {
    std::cout << "TEST B: Check Chemical Potential\n";
    std::cout << std::string(70, '-') << "\n";

    ThermoClass thermo;
    thermo.loadDatabase("CO.dat");
    thermo.setStandardUnits();

    // Test 001 conditions
    thermo.setTemperaturePressure(900.0, 1.0);
    thermo.setElementMass("C", 1.0);
    thermo.setElementMass("O", 1.0);

    int result = thermo.calculate();

    if (result != 0) {
        std::cout << "❌ Calculation failed: " << thermo.getErrorMessage() << "\n";
        return;
    }

    auto [muC, infoC] = thermo.getOutputChemPot("C");
    auto [muO, infoO] = thermo.getOutputChemPot("O");

    std::cout << "Test 001 (900K, 1 atm, C=1, O=1):\n";
    std::cout << "  μ_C = " << std::setprecision(6) << muC << " J/mol\n";
    std::cout << "  μ_O = " << muO << " J/mol\n";
    std::cout << "  Gibbs = " << thermo.getGibbsEnergy() << " J\n\n";

    // Try to get graphite moles
    auto [molesGraphite, infoGraphite] = thermo.getMolesPhase("C_Graphite(s)");
    std::cout << "  Graphite moles: " << molesGraphite << " mol\n";
    std::cout << "  Graphite query status: " << infoGraphite << "\n\n";

    // Check gas phase
    auto [molesGas, infoGas] = thermo.getMolesPhase("gas_ideal");
    std::cout << "  Gas moles: " << molesGas << " mol\n";
    std::cout << "  Gas query status: " << infoGas << "\n\n";

    std::cout << "Analysis:\n";
    if (muC > 0) {
        std::cout << "  ⚠️ μ_C > 0 indicates graphite SHOULD be stable!\n";
        std::cout << "     (Carbon chemical potential above standard state)\n";
    } else {
        std::cout << "  ✓ μ_C < 0 (graphite not thermodynamically favored)\n";
    }

    std::cout << "\n✓ Test B complete\n";
}

void testC_ForceGraphiteConstraint() {
    std::cout << "TEST C: Force Graphite Constraint\n";
    std::cout << std::string(70, '-') << "\n";

    ThermoClass thermo;
    thermo.loadDatabase("CO.dat");
    thermo.setStandardUnits();
    thermo.setTemperaturePressure(900.0, 1.0);
    thermo.setElementMass("C", 1.0);
    thermo.setElementMass("O", 1.0);

    std::cout << "Attempting to set constraint on C_Graphite(s) to 0.4 mol...\n";

    try {
        // Try the condensed phase constraint
        thermo.setCondPhaseConstraint("C_Graphite(s)", 0.4);
        std::cout << "✓ Constraint set successfully\n\n";

        int result = thermo.calculate();

        if (result != 0) {
            std::cout << "❌ Calculation failed: " << thermo.getErrorMessage() << "\n";
            std::cout << "   Error code: " << result << "\n";
        } else {
            std::cout << "✓ Calculation succeeded with constraint\n\n";

            auto [molesGraphite, info] = thermo.getMolesPhase("C_Graphite(s)");
            auto [molesGas, infoGas] = thermo.getMolesPhase("gas_ideal");

            std::cout << "Results:\n";
            std::cout << "  Graphite: " << molesGraphite << " mol\n";
            std::cout << "  Gas: " << molesGas << " mol\n";
            std::cout << "  Gibbs: " << thermo.getGibbsEnergy() << " J\n";

            if (std::abs(molesGraphite - 0.4) < 0.01) {
                std::cout << "\n✓ Constraint respected!\n";
            } else {
                std::cout << "\n⚠️ Constraint NOT respected (expected 0.4 mol)\n";
            }
        }
    } catch (const std::exception& e) {
        std::cout << "❌ Exception: " << e.what() << "\n";
    }

    std::cout << "\n✓ Test C complete\n";
}

void testD_PureGraphite() {
    std::cout << "TEST D: Pure Carbon System (No Oxygen)\n";
    std::cout << std::string(70, '-') << "\n";

    ThermoClass thermo;
    thermo.loadDatabase("CO.dat");
    thermo.setStandardUnits();
    thermo.setTemperaturePressure(900.0, 1.0);
    thermo.setElementMass("C", 1.0);
    // NO oxygen

    std::cout << "Running calculation with C=1 mol, O=0 mol...\n\n";

    int result = thermo.calculate();

    if (result != 0) {
        std::cout << "❌ Calculation failed: " << thermo.getErrorMessage() << "\n";
        std::cout << "   Error code: " << result << "\n";
    } else {
        std::cout << "✓ Calculation succeeded\n\n";

        auto [molesGraphite, infoG] = thermo.getMolesPhase("C_Graphite(s)");
        auto [molesGas, infoGas] = thermo.getMolesPhase("gas_ideal");

        std::cout << "Results:\n";
        std::cout << "  Graphite: " << molesGraphite << " mol\n";
        std::cout << "  Gas: " << molesGas << " mol\n";
        std::cout << "  Gibbs: " << thermo.getGibbsEnergy() << " J\n\n";

        std::cout << "Analysis:\n";
        if (molesGraphite > 0.9) {
            std::cout << "  ✓ Correct: Pure carbon → graphite\n";
        } else if (molesGas > 0.9) {
            std::cout << "  ❌ WRONG: Pure carbon → gas phase!\n";
            std::cout << "     (This is thermodynamically impossible)\n";
        } else {
            std::cout << "  ⚠️ Mixed phases (unexpected)\n";
        }
    }

    std::cout << "\n✓ Test D complete\n";
}

void testE_DetailedState() {
    std::cout << "TEST E: Detailed State Inspection\n";
    std::cout << std::string(70, '-') << "\n";

    ThermoClass thermo;
    thermo.loadDatabase("CO.dat");
    thermo.setStandardUnits();
    thermo.setTemperaturePressure(900.0, 1.0);
    thermo.setElementMass("C", 1.0);
    thermo.setElementMass("O", 1.0);
    thermo.calculate();

    const auto& ctx = thermo.getContext();
    const auto* state = ctx.thermo.get();

    std::cout << "System State:\n";
    std::cout << "  nElements: " << state->nElements << "\n";
    std::cout << "  nSolnPhasesSys: " << state->nSolnPhasesSys << "\n";
    std::cout << "  nConPhases: " << state->nConPhases << "\n";
    std::cout << "  nSpecies: " << state->nSpecies << "\n\n";

    std::cout << "Solution Phases in System:\n";
    for (int i = 0; i < state->nSolnPhasesSys; ++i) {
        std::cout << "  [" << i << "] " << state->cSolnPhaseName[i]
                  << " (moles: " << state->dMolesPhase[i] << ")\n";
    }

    std::cout << "\nPure Condensed Species:\n";
    std::cout << "  Total: " << state->nConPhases << "\n";
    for (int i = 0; i < std::min(state->nConPhases, 10); ++i) {
        int speciesIdx = state->nSpeciesPhase[0] + i;
        if (speciesIdx < state->nSpecies) {
            double moles = state->dMolFraction[speciesIdx];
            std::cout << "  [" << i << "] " << state->cSpeciesName[speciesIdx]
                      << " (moles: " << moles << ")\n";
        }
    }

    std::cout << "\nAssemblage:\n";
    std::cout << "  nConPhases in assemblage: " << state->nConPhases << "\n";
    for (int i = 0; i < state->nSolnPhasesSys; ++i) {
        int phaseID = state->iAssemblage[i];
        std::cout << "  Assemblage[" << i << "] = phase " << phaseID;
        if (phaseID > 0 && phaseID <= state->nSolnPhasesSys) {
            std::cout << " (" << state->cSolnPhaseName[phaseID-1] << ")";
        }
        std::cout << "\n";
    }

    std::cout << "\n✓ Test E complete\n";
}

int main() {
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║         GRAPHITE PHASE STABILITY DIAGNOSTIC TESTS                ║\n";
    std::cout << "║         Reference: PHASE1_DIAGNOSTICS.md                         ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n";

    printSeparator();
    testA_CheckGraphiteInDatabase();

    printSeparator();
    testB_CheckChemicalPotential();

    printSeparator();
    testC_ForceGraphiteConstraint();

    printSeparator();
    testD_PureGraphite();

    printSeparator();
    testE_DetailedState();

    printSeparator();
    std::cout << "All diagnostic tests complete!\n";
    std::cout << "Review results above to identify root cause.\n\n";

    return 0;
}
