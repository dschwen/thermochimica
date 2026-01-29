#pragma once

#include <Eigen/Dense>
#include <array>
#include "../util/Constants.hpp"

namespace Thermochimica {

/// Reinitialization state - converted from ModuleReinit.f90
/// Store previous computation state for warm restarts
struct ReinitState {
    // Previous assemblage
    Eigen::VectorXi iAssemblage_Old;

    // Elements used (fixed size 0:168)
    std::array<int, Constants::kMaxIsotopes> iElementsUsed_Old{};

    // Previous state vectors
    Eigen::VectorXd dChemicalPotential_Old;
    Eigen::VectorXd dMolesPhase_Old;
    Eigen::VectorXd dElementPotential_Old;
    Eigen::VectorXd dMolFraction_Old;

    // Constructor
    ReinitState() = default;

    /// Allocate arrays
    void allocate(int nElements, int nSpecies);

    /// Save current state for reinitialization
    void saveState(const Eigen::VectorXi& assemblage,
                   const Eigen::VectorXd& chemPot,
                   const Eigen::VectorXd& molesPhase,
                   const Eigen::VectorXd& elemPot,
                   const Eigen::VectorXd& molFrac);

    /// Check if reinit data is valid
    bool isValid() const;

    /// Clear saved state
    void clear();
};

} // namespace Thermochimica
