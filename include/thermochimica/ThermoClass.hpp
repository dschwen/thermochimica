/// @file ThermoClass.hpp
/// @brief Main Thermochimica class for object-oriented API
/// @details Class-based interface that replaces free-function API.
/// Owns all strategy instances and coordinates computation.

#pragma once

#include "thermochimica/ThermoContext.hpp"
#include <string>
#include <memory>
#include <vector>

namespace Thermochimica {

// Forward declarations
class ISolver;
class INewtonSolver;
class ILineSearch;
class IParser;
class IThermodynamicModel;
class PhaseAssemblageManager;
class ModelFactory;
struct ThermoState;
struct ThermoIO;
struct GEMState;
struct PhaseConstraints;

/// @brief Main Thermochimica class for thermodynamic calculations
/// @details Object-oriented interface using composition and strategy patterns.
/// Replaces the free-function API that takes ThermoContext& as first parameter.
///
/// Example usage:
/// @code
/// Thermochimica thermo;
/// thermo.loadDatabase("database.dat");
/// thermo.setTemperaturePressure(1000.0, 1.0);
/// thermo.setElementMass(6, 1.0);  // Carbon
/// thermo.calculate();
/// auto moles = thermo.getMolesPhase("Graphite");
/// @endcode
class ThermoClass {
public:
    /// @brief Constructor - initializes with default strategies
    ThermoClass();

    /// @brief Destructor
    ~ThermoClass();

    /// @brief Move constructor
    ThermoClass(ThermoClass&&) noexcept;

    /// @brief Move assignment
    ThermoClass& operator=(ThermoClass&&) noexcept;

    // Delete copy operations (ThermoContext is not copyable)
    ThermoClass(const ThermoClass&) = delete;
    ThermoClass& operator=(const ThermoClass&) = delete;

    // =========================================================================
    // Database Loading
    // =========================================================================

    /// @brief Set database filename
    /// @param filename Path to thermodynamic database file
    void setThermoFilename(const std::string& filename);

    /// @brief Parse database file
    /// @return Error code (0 = success)
    int parseCSDataFile();

    /// @brief Load database (combines setThermoFilename + parseCSDataFile)
    /// @param filename Path to database file
    /// @return Error code (0 = success)
    int loadDatabase(const std::string& filename);

    /// @brief Set custom parser
    /// @param parser Parser implementation
    void setParser(std::unique_ptr<IParser> parser);

    // =========================================================================
    // Input Configuration
    // =========================================================================

    /// @brief Set temperature in Kelvin
    /// @param temperature Temperature [K]
    void setTemperature(double temperature);

    /// @brief Set pressure in atm
    /// @param pressure Pressure [atm]
    void setPressure(double pressure);

    /// @brief Set temperature and pressure
    /// @param temperature Temperature [K]
    /// @param pressure Pressure [atm]
    void setTemperaturePressure(double temperature, double pressure);

    /// @brief Set element mass by atomic number
    /// @param atomicNumber Atomic number (e.g., 6 for carbon)
    /// @param mass Element mass [mass units]
    void setElementMass(int atomicNumber, double mass);

    /// @brief Set element mass by name
    /// @param elementName Element name (e.g., "C", "O")
    /// @param mass Element mass [mass units]
    void setElementMass(const std::string& elementName, double mass);

    /// @brief Set standard units (K, atm, moles)
    void setStandardUnits();

    /// @brief Set units to SI (K, Pa, moles)
    void setUnitsSI();

    /// @brief Set all units (temperature, pressure, mass)
    /// @param tempUnit Temperature unit (e.g., "K", "C", "F")
    /// @param pressUnit Pressure unit (e.g., "atm", "Pa", "bar")
    /// @param massUnit Mass unit (e.g., "moles", "grams", "kilograms")
    void setUnits(const std::string& tempUnit, const std::string& pressUnit, const std::string& massUnit);

    /// @brief Set temperature unit
    /// @param unit Temperature unit (e.g., "K", "C", "F")
    void setUnitTemperature(const std::string& unit);

    /// @brief Set pressure unit
    /// @param unit Pressure unit (e.g., "atm", "Pa", "bar")
    void setUnitPressure(const std::string& unit);

    /// @brief Set mass unit
    /// @param unit Mass unit (e.g., "moles", "grams", "kilograms")
    void setUnitMass(const std::string& unit);

    /// @brief Enable/disable fuzzy stoichiometry
    /// @param enable True to enable fuzzy stoichiometry
    void setFuzzyStoich(bool enable);

    /// @brief Set fuzzy stoichiometry magnitude
    /// @param magnitude Fuzzy stoichiometry magnitude threshold
    void setFuzzyMagnitude(double magnitude);

    /// @brief Set print results mode
    /// @param mode Print mode (0 = none, 1 = summary, 2 = detailed)
    void setPrintResultsMode(int mode);

    /// @brief Enable/disable JSON output writing
    /// @param enable True to enable JSON output
    void setWriteJSON(bool enable);

    /// @brief Enable heat capacity, entropy, and enthalpy calculations
    /// @param enable True to enable HCP/S/H calculations
    void setHeatCapacityEntropyEnthalpy(bool enable);

    // =========================================================================
    // Solver Configuration (Advanced)
    // =========================================================================

    /// @brief Set custom solver strategy
    /// @param solver Solver implementation
    void setSolver(std::unique_ptr<ISolver> solver);

    /// @brief Set custom Newton solver
    /// @param newton Newton solver implementation
    void setNewtonSolver(std::unique_ptr<INewtonSolver> newton);

    /// @brief Set custom line search
    /// @param lineSearch Line search implementation
    void setLineSearch(std::unique_ptr<ILineSearch> lineSearch);

    // =========================================================================
    // Phase Constraints
    // =========================================================================

    /// @brief Set solution phase constraint
    /// @param phaseName Phase name
    /// @param targetFraction Target element fraction [0-1]
    void setSolnPhaseConstraint(const std::string& phaseName, double targetFraction);

    /// @brief Set condensed phase constraint
    /// @param speciesName Species name
    /// @param targetFraction Target element fraction [0-1]
    void setCondPhaseConstraint(const std::string& speciesName, double targetFraction);

    /// @brief Remove phase constraint
    /// @param phaseName Phase name
    void removePhaseConstraint(const std::string& phaseName);

    /// @brief Clear all phase constraints
    void clearPhaseConstraints();

    /// @brief Set phase constraint by index
    /// @param phaseIndex Phase index in system
    /// @param isSolutionPhase True if solution phase, false if pure phase
    /// @param targetFraction Target element fraction [0-1]
    void setPhaseConstraint(int phaseIndex, bool isSolutionPhase, double targetFraction);

    /// @brief Get element fraction in phase
    /// @param phaseName Phase name
    /// @return Pair of (fraction, error code)
    std::pair<double, int> getPhaseElementFraction(const std::string& phaseName) const;

    /// @brief Check if phase constraints are satisfied
    /// @return True if all constraints satisfied
    bool arePhaseConstraintsSatisfied() const;

    /// @brief Set constraint tolerance
    /// @param tolerance Convergence tolerance for constrained solver
    void setConstraintTolerance(double tolerance);

    /// @brief Set constraint penalty parameter
    /// @param rho Penalty parameter (augmented Lagrangian)
    void setConstraintPenaltyParameter(double rho);

    /// @brief Set maximum constraint outer iterations
    /// @param maxIter Maximum outer iterations for constrained solver
    void setConstraintMaxOuterIterations(int maxIter);

    // =========================================================================
    // Main Computation
    // =========================================================================

    /// @brief Run thermodynamic calculation
    /// @details Replaces thermochimica(ctx) free function
    /// @return Error code (0 = success)
    int calculate();

    // =========================================================================
    // Granular Control (Advanced)
    // =========================================================================

    /// @brief Initialize calculation
    void initialize();

    /// @brief Check system validity
    void checkSystem();

    /// @brief Compute thermodynamic data
    void computeThermoData();

    /// @brief Setup for iteration
    void setup();

    /// @brief Run solver (after setup)
    /// @return Error code (0 = success)
    int solve();

    // =========================================================================
    // Output Retrieval
    // =========================================================================

    /// @brief Get chemical potential
    /// @param elementName Element name
    /// @return Pair of (chemical potential, error code)
    std::pair<double, int> getOutputChemPot(const std::string& elementName) const;

    /// @brief Get phase moles
    /// @param phaseName Phase name
    /// @return Pair of (moles, error code)
    std::pair<double, int> getMolesPhase(const std::string& phaseName) const;

    /// @brief Get system Gibbs energy
    /// @return Gibbs energy [J]
    double getGibbsEnergy() const;

    /// @brief Get heat capacity
    /// @return Heat capacity [J/K]
    double getHeatCapacity() const;

    /// @brief Get entropy
    /// @return Entropy [J/K]
    double getEntropy() const;

    /// @brief Get enthalpy
    /// @return Enthalpy [J]
    double getEnthalpy() const;

    /// @brief Get species mole fraction and moles in solution phase
    /// @param phaseName Solution phase name
    /// @param speciesName Species name
    /// @return Tuple of (mole fraction, moles, error code)
    std::tuple<double, double, int> getOutputSolnSpecies(const std::string& phaseName,
                                                          const std::string& speciesName) const;

    /// @brief Get species mole fraction and moles (alternative method)
    /// @param phaseName Phase name
    /// @param speciesName Species name
    /// @return Tuple of (mole fraction, moles, error code)
    std::tuple<double, double, int> getOutputMolSpecies(const std::string& phaseName,
                                                         const std::string& speciesName) const;

    /// @brief Get moles of solution phase
    /// @param phaseName Solution phase name
    /// @return Pair of (moles, error code)
    std::pair<double, int> getSolnPhaseMol(const std::string& phaseName) const;

    /// @brief Get moles of pure condensed phase
    /// @param phaseName Pure condensed phase name
    /// @return Pair of (moles, error code)
    std::pair<double, int> getPureConPhaseMol(const std::string& phaseName) const;

    /// @brief Get moles of element in phase
    /// @param elementName Element name
    /// @param phaseName Phase name
    /// @return Pair of (moles, error code)
    std::pair<double, int> getElementMolesInPhase(const std::string& elementName,
                                                   const std::string& phaseName) const;

    /// @brief Get phase index by name
    /// @param phaseName Phase name
    /// @return Phase index (-1 if not found)
    int getPhaseIndex(const std::string& phaseName) const;

    /// @brief Get site fraction for sublattice model
    /// @param phaseName Phase name
    /// @param sublattice Sublattice index
    /// @param constituentName Constituent name
    /// @return Pair of (site fraction, error code)
    std::pair<double, int> getOutputSiteFraction(const std::string& phaseName,
                                                  int sublattice,
                                                  const std::string& constituentName) const;

    /// @brief Check if phase is gas phase
    /// @param phaseName Phase name
    /// @return True if gas phase
    bool isPhaseGas(const std::string& phaseName) const;

    /// @brief Check if phase is MQM (Modified Quasichemical Model)
    /// @param phaseName Phase name
    /// @return True if MQM phase
    bool isPhaseMQM(const std::string& phaseName) const;

    /// @brief Get element chemical potential by index
    /// @param elementIndex Element index
    /// @return Pair of (chemical potential, error code)
    std::pair<double, int> getElementChemicalPotential(int elementIndex) const;

    /// @brief Get all element chemical potentials
    /// @return Vector of chemical potentials
    std::vector<double> getAllElementChemicalPotentials() const;

    /// @brief Get Gibbs energy derivative for phase
    /// @param phaseName Phase name
    /// @return Pair of (derivative, error code)
    std::pair<double, int> getGibbsEnergyDerivative(const std::string& phaseName) const;

    /// @brief Print results to stdout
    void printResults() const;

    /// @brief Print detailed results to stdout
    void printResultsDetailed() const;

    /// @brief Write results to JSON file
    /// @param append If true, append to existing file
    void writeJSON(bool append = false);

    /// @brief Post-process calculation results
    void postProcess();

    // =========================================================================
    // Status
    // =========================================================================

    /// @brief Get error/info code
    /// @return Error code (0 = success)
    int getInfoCode() const;

    /// @brief Check if calculation succeeded
    /// @return true if successful
    bool isSuccess() const;

    /// @brief Get error message
    /// @return Error description
    std::string getErrorMessage() const;

    // =========================================================================
    // Reset
    // =========================================================================

    /// @brief Reset for new calculation (keeps database)
    void reset();

    /// @brief Reset everything including database
    void resetAll();

    // =========================================================================
    // Database Queries
    // =========================================================================

    /// @brief Get number of elements in database
    /// @return Number of elements
    int getNumberElementsDatabase() const;

    /// @brief Get element name at index
    /// @param index Element index (0-based)
    /// @return Element name
    std::string getElementAtIndex(int index) const;

    /// @brief Get number of phases in system
    /// @return Pair of (solution phases, pure condensed phases)
    std::pair<int, int> getNumberPhasesSystem() const;

    /// @brief Get phase name at index
    /// @param index Phase index
    /// @return Phase name
    std::string getPhaseNameAtIndex(int index) const;

    /// @brief Get number of species in system
    /// @return Number of species
    int getNumberSpeciesSystem() const;

    /// @brief Get species name at index
    /// @param index Species index
    /// @return Species name
    std::string getSpeciesAtIndex(int index) const;

    // =========================================================================
    // Reinitialization
    // =========================================================================

    /// @brief Save current state for reinitialization
    void saveReinitData();

    /// @brief Set reinit requested flag
    /// @param requested True to request reinitialization
    void setReinitRequested(bool requested);

    /// @brief Check if reinit data is available
    /// @return True if reinit data available
    bool isReinitDataAvailable() const;

    /// @brief Get reinitialization data
    /// @return Tuple of (assemblage, species moles, mole fractions, chemical potentials, element potentials)
    std::tuple<std::vector<int>, std::vector<double>, std::vector<double>,
               std::vector<double>, std::vector<double>> getReinitData() const;

    /// @brief Set reinitialization data
    /// @param assemblage Assemblage vector
    /// @param speciesMoles Species moles vector
    /// @param moleFractions Mole fractions vector
    /// @param chemicalPotentials Chemical potentials vector
    /// @param elementPotentials Element potentials vector
    void setReinitData(const std::vector<int>& assemblage,
                       const std::vector<double>& speciesMoles,
                       const std::vector<double>& moleFractions,
                       const std::vector<double>& chemicalPotentials,
                       const std::vector<double>& elementPotentials);

    // =========================================================================
    // Utility Functions (Static)
    // =========================================================================

    /// @brief Get atomic number from element symbol
    /// @param symbol Element symbol (e.g., "C", "O")
    /// @return Atomic number (-1 if not found)
    static int getAtomicNumber(const std::string& symbol);

    /// @brief Get element symbol from atomic number
    /// @param atomicNumber Atomic number
    /// @return Element symbol (empty if invalid)
    static std::string getElementSymbol(int atomicNumber);

    // =========================================================================
    // Legacy Compatibility
    // =========================================================================

    /// @brief Get underlying context (for legacy code)
    /// @return Reference to ThermoContext
    ThermoContext& getContext() { return context_; }

    /// @brief Get underlying context (const)
    /// @return Const reference to ThermoContext
    const ThermoContext& getContext() const { return context_; }

private:
    // Core state (owns all data)
    ThermoContext context_;

    // Convenience pointers (non-owning)
    ThermoState* state_;
    ThermoIO* io_;
    GEMState* gemState_;
    PhaseConstraints* phaseConstraints_;

    // Strategy components (owned)
    std::unique_ptr<IParser> parser_;
    std::unique_ptr<ISolver> solver_;
    std::unique_ptr<INewtonSolver> newtonSolver_;
    std::unique_ptr<ILineSearch> lineSearch_;

    // Composition components (owned)
    std::unique_ptr<PhaseAssemblageManager> phaseManager_;
    std::unique_ptr<ModelFactory> modelFactory_;

    /// @brief Initialize default strategies
    void initializeDefaultStrategies();

    /// @brief Validate configuration before calculation
    void validateConfiguration() const;
};

} // namespace Thermochimica
