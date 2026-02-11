/// @file CouplingAPI.h
/// @brief C API for Thermochimica - for use from C code or other languages via FFI
/// @details Provides a simple C interface to Thermochimica thermodynamic calculations.
/// Uses a global singleton instance for backward compatibility.

#ifndef THERMOCHIMICA_COUPLING_API_H
#define THERMOCHIMICA_COUPLING_API_H

#ifdef __cplusplus
extern "C" {
#endif

/// @brief Initialize the Thermochimica instance
/// @details Creates the global ThermoClass instance if not already created.
/// Most functions call this automatically, so explicit calls are usually unnecessary.
void TC_init();

/// @brief Clean up and destroy the Thermochimica instance
/// @details Frees all allocated memory. Call this before program exit.
void TC_cleanup();

/// @brief Load thermodynamic database from file
/// @param filename Path to database file (e.g., "CO.dat")
/// @details This is the recommended function to load a database.
/// Combines setting the filename and parsing in one step.
void TC_loadDatabase(const char* filename);

/// @brief Set thermodynamic database filename (DEPRECATED)
/// @param filename Path to database file
/// @deprecated Use TC_loadDatabase() instead. This function immediately loads
/// the database, despite its name suggesting it only sets a filename.
void TC_setThermoFilename(const char* filename);

/// @brief Parse database file (DEPRECATED - NO-OP)
/// @deprecated This function does nothing. Use TC_loadDatabase() instead,
/// which loads and parses the database in one step.
void TC_parseDataFile();

/// @brief Set temperature and pressure for calculation
/// @param T Temperature [K]
/// @param P Pressure [atm]
void TC_setTemperaturePressure(double T, double P);

/// @brief Set element mass/amount for calculation
/// @param atomicNum Atomic number (e.g., 6 for Carbon, 8 for Oxygen)
/// @param mass Element mass/amount [moles by default]
void TC_setElementMass(int atomicNum, double mass);

/// @brief Run thermodynamic calculation
/// @details Performs Gibbs energy minimization to find equilibrium phase assemblage.
/// Check TC_getInfo() for success/error status after calling.
void TC_solve();

/// @brief Get calculation status/error code
/// @return 0 on success, non-zero error code on failure, -1 if not initialized
int TC_getInfo();

/// @brief Get system Gibbs energy from last calculation
/// @return Gibbs energy [J], or 0.0 if not initialized
double TC_getGibbsEnergy();

/// @brief Reset for new calculation (keeps database loaded)
/// @details Clears calculation results but preserves loaded database.
void TC_reset();

/// @brief Reset everything including database
/// @details Clears all data including loaded database.
void TC_resetAll();

#ifdef __cplusplus
}
#endif

#endif // THERMOCHIMICA_COUPLING_API_H
