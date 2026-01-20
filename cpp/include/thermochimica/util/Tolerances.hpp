#pragma once

#include <array>
#include "Constants.hpp"

namespace Thermochimica {

// Tolerance indices (matching Fortran indexing for clarity)
enum ToleranceIndex {
    kTolMassBalance = 0,       // Mass balance residual
    kTolChemPotential = 1,     // Chemical potential residual
    kTolGibbsEnergy = 2,       // Gibbs energy tolerance
    kTolMoleFraction = 3,      // Mole fraction tolerance
    kTolDrivingForce = 4,      // Driving force tolerance
    kTolElementPotential = 5,  // Element potential tolerance
    kTolConvergence = 6,       // Convergence tolerance
    kTolPhaseChange = 7,       // Phase change tolerance
    kTolSiteFraction = 8,      // Site fraction tolerance
    kTolMinMoles = 9,          // Minimum moles threshold
    kTolPhaseMoles = 10,       // Phase moles tolerance
    kTolFunctionNorm = 11,     // Function norm tolerance
    kTolRelative = 12,         // Relative tolerance
    kTolAbsolute = 13,         // Absolute tolerance
    kTolEuclidean = 14         // Euclidean norm tolerance
};

struct Tolerances {
    std::array<double, Constants::kNumTolerances> values;

    Tolerances() {
        initDefaults();
    }

    void initDefaults() {
        // Default tolerance values (from InitThermo.f90)
        values[kTolMassBalance] = 1.0e-5;
        values[kTolChemPotential] = 1.0e-6;
        values[kTolGibbsEnergy] = 1.0e-8;
        values[kTolMoleFraction] = 1.0e-12;
        values[kTolDrivingForce] = 1.0e-4;
        values[kTolElementPotential] = 1.0e-6;
        values[kTolConvergence] = 1.0e-6;
        values[kTolPhaseChange] = 1.0e-3;
        values[kTolSiteFraction] = 1.0e-6;
        values[kTolMinMoles] = 1.0e-20;
        values[kTolPhaseMoles] = 1.0e-9;
        values[kTolFunctionNorm] = 1.0e-6;
        values[kTolRelative] = 1.0e-8;
        values[kTolAbsolute] = 1.0e-10;
        values[kTolEuclidean] = 1.0e-2;
    }

    double& operator[](int index) { return values[index]; }
    const double& operator[](int index) const { return values[index]; }
};

} // namespace Thermochimica
