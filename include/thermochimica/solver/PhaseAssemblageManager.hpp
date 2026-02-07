/// @file PhaseAssemblageManager.hpp
/// @brief Manager for phase assemblage operations
/// @details Encapsulates phase add/remove logic and driving force computation

#pragma once

namespace Thermochimica {

// Forward declarations
class ThermoState;
class GEMState;
class ThermoIO;

/// @brief Manages phase assemblage changes during GEM iterations
/// @details Handles adding/removing phases based on driving forces,
/// computes driving forces for unstable phases, and maintains assemblage history
class PhaseAssemblageManager {
public:
    /// @brief Constructor
    /// @param state Reference to thermodynamic state
    /// @param gemState Reference to GEM solver state
    /// @param io Reference to ThermoIO for temperature
    PhaseAssemblageManager(ThermoState& state, GEMState& gemState, ThermoIO& io);

    /// @brief Check and adjust phase assemblage
    /// @details Main entry point - computes driving forces and adds/removes phases
    void check();

    /// @brief Add a solution phase to the assemblage
    /// @param phaseIndex System phase index (0-based)
    /// @return true if phase was added successfully
    bool addSolnPhase(int phaseIndex);

    /// @brief Remove a solution phase from the assemblage
    /// @param phaseIndex Assemblage index for the phase
    /// @return true if phase was removed successfully
    bool removeSolnPhase(int phaseIndex);

    /// @brief Add a pure condensed phase to the assemblage
    /// @param speciesIndex Species index (0-based)
    /// @return true if phase was added successfully
    bool addPureConPhase(int speciesIndex);

    /// @brief Remove a pure condensed phase from the assemblage
    /// @param speciesIndex Assemblage index for the phase
    /// @return true if phase was removed successfully
    bool removePureConPhase(int speciesIndex);

    /// @brief Revert to previous assemblage
    /// @details Restores saved assemblage state
    void revert();

private:
    /// @brief Check if a species is feasible (contains only active elements)
    /// @param speciesIdx Species index
    /// @return true if species contains only elements with non-zero input mass
    bool speciesIsFeasible(int speciesIdx) const;

    /// @brief Compute driving forces for all unstable solution phases
    /// @details For each solution phase not in assemblage, estimates optimal
    /// mole fractions and computes driving force
    void computeSolutionPhaseDrivingForces();

    /// @brief Compute driving forces for all unstable pure condensed phases
    /// @details For each pure phase not in assemblage, computes driving force
    void computePureConPhaseDrivingForces();

    /// @brief Find phase with maximum driving force
    /// @param isBestSoln Output: true if best is solution phase, false if pure con
    /// @param bestIndex Output: phase/species index
    /// @return Maximum driving force value
    double findMaxDrivingForce(bool& isBestSoln, int& bestIndex) const;

    /// @brief Save current assemblage for potential revert
    void saveCurrentAssemblage();

    ThermoState& state_;   ///< Reference to thermodynamic state
    GEMState& gemState_;   ///< Reference to GEM solver state
    ThermoIO& io_;         ///< Reference to ThermoIO
};

} // namespace Thermochimica
