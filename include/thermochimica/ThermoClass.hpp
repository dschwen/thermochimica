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

    /// @brief Print results to stdout
    void printResults() const;

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
