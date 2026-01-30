#pragma once

#include "ThermoContext.hpp"
#include "util/Constants.hpp"
#include "util/ErrorCodes.hpp"
#include "util/Tolerances.hpp"

#include <string>
#include <tuple>
#include <vector>

namespace Thermochimica {

// ============================================================================
// Main Solver Functions
// ============================================================================

/// Main entry point - performs complete equilibrium calculation
/// Equivalent to the Fortran Thermochimica() subroutine
void thermochimica(ThermoContext& ctx);

/// Initialize the solver
void init(ThermoContext& ctx);

/// Check system validity
void checkSystem(ThermoContext& ctx);

/// Compute thermodynamic data
void compThermoData(ThermoContext& ctx);

/// Setup phase assemblage
void setup(ThermoContext& ctx);

/// Run the GEM solver
void solve(ThermoContext& ctx);

// ============================================================================
// Parser Functions
// ============================================================================

/// Parse a ChemSage data file
void parseCSDataFile(ThermoContext& ctx, const std::string& filename);

/// Parse using filename set in context
void parseCSDataFile(ThermoContext& ctx);

// ============================================================================
// Input Setting Functions
// ============================================================================

/// Set the thermodynamic data filename
void setThermoFilename(ThermoContext& ctx, const std::string& filename);

/// Set temperature and pressure
void setTemperaturePressure(ThermoContext& ctx, double temperature, double pressure);

/// Set temperature
void setTemperature(ThermoContext& ctx, double temperature);

/// Set pressure
void setPressure(ThermoContext& ctx, double pressure);

/// Set element mass by atomic number (1-118)
void setElementMass(ThermoContext& ctx, int atomicNumber, double mass);

/// Set element mass by name
void setElementMass(ThermoContext& ctx, const std::string& elementName, double mass);

/// Set units (temperature, pressure, mass)
void setUnits(ThermoContext& ctx,
              const std::string& tempUnit,
              const std::string& pressUnit,
              const std::string& massUnit);

/// Set temperature unit ('K', 'C', 'F', 'R')
void setUnitTemperature(ThermoContext& ctx, const std::string& unit);

/// Set pressure unit ('atm', 'psi', 'bar', 'Pa', 'kPa')
void setUnitPressure(ThermoContext& ctx, const std::string& unit);

/// Set mass unit ('moles', 'grams', 'kilograms', etc.)
void setUnitMass(ThermoContext& ctx, const std::string& unit);

/// Set standard units (K, atm, moles)
void setStandardUnits(ThermoContext& ctx);

/// Enable/disable fuzzy stoichiometry
void setFuzzyStoich(ThermoContext& ctx, bool enable);

/// Set fuzzy magnitude
void setFuzzyMagnitude(ThermoContext& ctx, double magnitude);

/// Set print results mode (0=none, 1=basic, 2=detailed)
void setPrintResultsMode(ThermoContext& ctx, int mode);

/// Enable/disable JSON output
void setWriteJSON(ThermoContext& ctx, bool enable);

/// Enable/disable heat capacity, entropy, and enthalpy calculation
void setHeatCapacityEntropyEnthalpy(ThermoContext& ctx, bool enable);

// ============================================================================
// Phase Constraint Functions
// ============================================================================

/// Set phase fraction constraint for a solution phase by name
/// @param ctx The thermochimica context
/// @param phaseName Name of the solution phase
/// @param targetFraction Target element-level phase fraction [0, 1]
void setSolnPhaseConstraint(ThermoContext& ctx,
                            const std::string& phaseName,
                            double targetFraction);

/// Set phase fraction constraint for a pure condensed phase by species name
/// @param ctx The thermochimica context
/// @param speciesName Name of the pure condensed phase species
/// @param targetFraction Target element-level phase fraction [0, 1]
void setCondPhaseConstraint(ThermoContext& ctx,
                            const std::string& speciesName,
                            double targetFraction);

/// Set phase fraction constraint by phase index
/// @param ctx The thermochimica context
/// @param phaseIndex Index of the phase (0-based)
/// @param isSolutionPhase True for solution phase, false for pure condensed
/// @param targetFraction Target element-level phase fraction [0, 1]
void setPhaseConstraint(ThermoContext& ctx,
                        int phaseIndex,
                        bool isSolutionPhase,
                        double targetFraction);

/// Remove phase fraction constraint
/// @param ctx The thermochimica context
/// @param phaseName Name of the phase to remove constraint from
void removePhaseConstraint(ThermoContext& ctx,
                           const std::string& phaseName);

/// Remove all phase constraints
/// @param ctx The thermochimica context
void clearPhaseConstraints(ThermoContext& ctx);

/// Get current element-level phase fraction
/// The element phase fraction is: (sum of element moles in phase) / (total element moles)
/// @param ctx The thermochimica context
/// @param phaseName Name of the phase
/// @return (fraction, info) where info=0 on success
std::pair<double, int> getPhaseElementFraction(const ThermoContext& ctx,
                                               const std::string& phaseName);

/// Check if all phase constraints are satisfied within tolerance
/// @param ctx The thermochimica context
/// @return True if all constraints are satisfied
bool arePhaseConstraintsSatisfied(const ThermoContext& ctx);

/// Get all element chemical potentials (dimensionless, divided by RT)
/// @param ctx The thermochimica context
/// @return Vector of element chemical potentials indexed by element order in database
std::vector<double> getAllElementChemicalPotentials(const ThermoContext& ctx);

/// Get element chemical potential by index (dimensionless, divided by RT)
/// @param ctx The thermochimica context
/// @param elementIndex Element index (0-based)
/// @return (chemPot, info) where info=0 on success
std::pair<double, int> getElementChemicalPotential(const ThermoContext& ctx, int elementIndex);

/// Get derivative of Gibbs energy with respect to phase fraction (dG/df)
/// This is the driving force for phase fraction change at constant T, P, composition.
/// For constrained equilibrium, this equals -λ (negative Lagrange multiplier).
/// Units: J (same as Gibbs energy)
/// @param ctx The thermochimica context
/// @param phaseName Name of the phase
/// @return (dGdf, info) where info=0 on success
std::pair<double, int> getGibbsEnergyDerivative(const ThermoContext& ctx,
                                                 const std::string& phaseName);

/// Set constraint tolerance
/// @param ctx The thermochimica context
/// @param tolerance Tolerance for constraint satisfaction (default 1e-4)
void setConstraintTolerance(ThermoContext& ctx, double tolerance);

/// Set penalty parameter for augmented Lagrangian method
/// @param ctx The thermochimica context
/// @param rho Initial penalty parameter (default 1.0)
void setConstraintPenaltyParameter(ThermoContext& ctx, double rho);

/// Set maximum outer iterations for augmented Lagrangian loop
/// @param ctx The thermochimica context
/// @param maxIter Maximum number of outer iterations (default 20)
void setConstraintMaxOuterIterations(ThermoContext& ctx, int maxIter);

// ============================================================================
// Output Retrieval Functions
// ============================================================================

/// Get element chemical potential by name
/// Returns (chemPot, info) where info=0 on success
std::pair<double, int> getOutputChemPot(const ThermoContext& ctx,
                                        const std::string& elementName);

/// Get solution species output (mole fraction, chemical potential)
/// Returns (moleFrac, chemPot, info)
std::tuple<double, double, int> getOutputSolnSpecies(
    const ThermoContext& ctx,
    const std::string& phaseName,
    const std::string& speciesName);

/// Get species moles and mole fraction
/// Returns (moles, moleFrac, info)
std::tuple<double, double, int> getOutputMolSpecies(
    const ThermoContext& ctx,
    const std::string& phaseName,
    const std::string& speciesName);

/// Get moles of a phase by name
/// Returns (moles, info)
std::pair<double, int> getMolesPhase(const ThermoContext& ctx,
                                     const std::string& phaseName);

/// Get solution phase moles
/// Returns (moles, info)
std::pair<double, int> getSolnPhaseMol(const ThermoContext& ctx,
                                       const std::string& phaseName);

/// Get pure condensed phase moles
/// Returns (moles, info)
std::pair<double, int> getPureConPhaseMol(const ThermoContext& ctx,
                                          const std::string& phaseName);

/// Get element moles in a phase
/// Returns (moles, info)
std::pair<double, int> getElementMolesInPhase(
    const ThermoContext& ctx,
    const std::string& elementName,
    const std::string& phaseName);

/// Get phase index by name (-1 if not found)
int getPhaseIndex(const ThermoContext& ctx, const std::string& phaseName);

/// Get site fraction for sublattice phases
/// Returns (siteFrac, info)
std::pair<double, int> getOutputSiteFraction(
    const ThermoContext& ctx,
    const std::string& phaseName,
    int sublattice,
    const std::string& constituentName);

/// Check if phase is a gas phase
bool isPhaseGas(const ThermoContext& ctx, const std::string& phaseName);

/// Check if phase is MQM (modified quasichemical)
bool isPhaseMQM(const ThermoContext& ctx, const std::string& phaseName);

/// Get system Gibbs energy
double getGibbsEnergy(const ThermoContext& ctx);

/// Get heat capacity (requires lHeatCapacityEntropyEnthalpy = true)
double getHeatCapacity(const ThermoContext& ctx);

/// Get entropy
double getEntropy(const ThermoContext& ctx);

/// Get enthalpy
double getEnthalpy(const ThermoContext& ctx);

// ============================================================================
// Database Query Functions
// ============================================================================

/// Get number of elements in loaded database
int getNumberElementsDatabase(const ThermoContext& ctx);

/// Get element name at index (0-based)
std::string getElementAtIndex(const ThermoContext& ctx, int index);

/// Get number of phases (solution, pure condensed)
std::pair<int, int> getNumberPhasesSystem(const ThermoContext& ctx);

/// Get phase name at index (0-based)
std::string getPhaseNameAtIndex(const ThermoContext& ctx, int index);

/// Get number of species in system
int getNumberSpeciesSystem(const ThermoContext& ctx);

/// Get species name at index (0-based)
std::string getSpeciesAtIndex(const ThermoContext& ctx, int index);

// ============================================================================
// Reinitialization Functions
// ============================================================================

/// Save current state for reinitialization
void saveReinitData(ThermoContext& ctx);

/// Request reinitialization from saved data
void setReinitRequested(ThermoContext& ctx, bool requested);

/// Check if reinit data is available
bool isReinitDataAvailable(const ThermoContext& ctx);

/// Get reinitialization data
/// Returns vectors of (assemblage, chemPot, molesPhase, elemPot, molFrac)
std::tuple<std::vector<int>, std::vector<double>, std::vector<double>,
           std::vector<double>, std::vector<double>>
getReinitData(const ThermoContext& ctx);

/// Set reinitialization data
void setReinitData(ThermoContext& ctx,
                   const std::vector<int>& assemblage,
                   const std::vector<double>& chemPot,
                   const std::vector<double>& molesPhase,
                   const std::vector<double>& elemPot,
                   const std::vector<double>& molFrac);

// ============================================================================
// Reset Functions
// ============================================================================

/// Reset for new calculation (keeps database loaded)
void resetThermo(ThermoContext& ctx);

/// Full reset including database
void resetThermoAll(ThermoContext& ctx);

// ============================================================================
// Output Functions
// ============================================================================

/// Post-process results (compute Gibbs energy, prepare output arrays)
void postProcess(ThermoContext& ctx);

/// Compute heat capacity, entropy, and enthalpy
/// Requires io.lHeatCapacityEntropyEnthalpy = true
void computeHeatCapacity(ThermoContext& ctx);

/// Print results to stdout
void printResults(ThermoContext& ctx);

/// Print detailed results to stdout
void printResultsDetailed(ThermoContext& ctx);

/// Write JSON output
void writeJSON(ThermoContext& ctx, bool append = false);

// ============================================================================
// Utility Functions
// ============================================================================

/// Get atomic number from element symbol
int getAtomicNumber(const std::string& symbol);

/// Get element symbol from atomic number
std::string getElementSymbol(int atomicNumber);

/// Get error message for error code
const char* getErrorMessage(int errorCode);

} // namespace Thermochimica
