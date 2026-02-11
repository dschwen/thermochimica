/// Basic usage example for Thermochimica C++
/// Demonstrates the modern ThermoClass object-oriented API

#include <thermochimica/ThermoClass.hpp>
#include <iostream>

int main() {
    // Create a ThermoClass instance (RAII - automatic cleanup)
    Thermochimica::ThermoClass thermo;

    // Set units
    thermo.setStandardUnits();  // K, atm, moles

    // Load database (combines setThermoFilename + parseCSDataFile)
    int result = thermo.loadDatabase("CO.dat");

    if (result != 0) {
        std::cerr << "Error loading database (code " << result << "): "
                  << thermo.getErrorMessage()
                  << std::endl;
        return 1;
    }

    std::cout << "Database loaded successfully.\n";
    std::cout << "Elements: " << thermo.getNumberElementsDatabase() << "\n";
    std::cout << "Species: " << thermo.getNumberSpeciesSystem() << "\n";

    // Set temperature and pressure
    thermo.setTemperaturePressure(1000.0, 1.0);  // 1000 K, 1 atm

    // Set composition (by atomic number)
    thermo.setElementMass(6, 1.0);   // 1 mole of Carbon
    thermo.setElementMass(8, 2.0);   // 2 moles of Oxygen

    // Run the calculation
    result = thermo.calculate();

    // Check result
    if (thermo.isSuccess()) {
        std::cout << "\nCalculation successful!\n";
        std::cout << "System Gibbs Energy: "
                  << thermo.getGibbsEnergy()
                  << " J\n";

        // Print results
        thermo.printResults();
    } else {
        std::cerr << "Calculation failed (code " << result << "): "
                  << thermo.getErrorMessage()
                  << std::endl;

        // Debug info (accessing internal context for diagnostics)
        auto& ctx = thermo.getContext();
        std::cerr << "\nDebug info:\n";
        std::cerr << "  nSpecies in thermo state: " << ctx.thermo->nSpecies << "\n";
        std::cerr << "  nElements in thermo state: " << ctx.thermo->nElements << "\n";
        std::cerr << "  nSolnPhasesSys: " << ctx.thermo->nSolnPhasesSys << "\n";
        std::cerr << "  nConPhases: " << ctx.thermo->nConPhases << "\n";
        std::cerr << "  Temperature: " << ctx.io->dTemperature << " K\n";
        std::cerr << "  Pressure: " << ctx.io->dPressure << " atm\n";

        std::cerr << "\n  Element masses (from IO state):\n";
        for (int i = 0; i < 20; ++i) {
            if (ctx.io->dElementMass[i] > 0.0) {
                std::cerr << "    Element " << i << ": " << ctx.io->dElementMass[i] << " moles\n";
            }
        }

        std::cerr << "\n  Element moles (in thermo state):\n";
        for (int j = 0; j < ctx.thermo->nElements; ++j) {
            std::cerr << "    " << ctx.thermo->cElementName[j] << " (idx " << j << "): "
                      << ctx.thermo->dMolesElement(j) << " moles\n";
        }

        std::cerr << "\n  Phase assemblage (iAssemblage):\n";
        for (int i = 0; i < std::min(5, static_cast<int>(ctx.thermo->iAssemblage.size())); ++i) {
            std::cerr << "    Slot " << i << ": species " << ctx.thermo->iAssemblage(i) << "\n";
        }

        std::cerr << "\n  First 5 species info:\n";
        for (int i = 0; i < std::min(5, ctx.thermo->nSpecies); ++i) {
            std::cerr << "    Species " << i << " (" << ctx.thermo->cSpeciesName[i] << "): ";
            std::cerr << "G=" << ctx.thermo->dStdGibbsEnergy(i) << " J/mol, ";
            std::cerr << "mu/RT=" << ctx.thermo->dChemicalPotential(i) << "\n";
        }

        std::cerr << "\n  Stoichiometry matrix (first 5 species):\n";
        for (int i = 0; i < std::min(5, ctx.thermo->nSpecies); ++i) {
            std::cerr << "    " << ctx.thermo->cSpeciesName[i] << ": ";
            for (int j = 0; j < ctx.thermo->nElements; ++j) {
                std::cerr << ctx.thermo->dStoichSpecies(i, j) << " ";
            }
            std::cerr << "\n";
        }

        return 1;
    }

    return 0;
}
