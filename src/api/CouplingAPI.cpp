// C-compatible API for external coupling
// Provides C interface functions for use from C code or other languages via FFI

#include "thermochimica/ThermoClass.hpp"
#include "thermochimica/util/ErrorCodes.hpp"
#include <exception>

extern "C" {

// Global ThermoClass instance for C API (singleton pattern for backward compatibility)
static Thermochimica::ThermoClass* g_thermo = nullptr;

// Exception-safe wrapper macro for C API functions
#define TC_CATCH_EXCEPTIONS(error_code) \
    catch (const std::exception& e) { \
        if (g_thermo) { \
            g_thermo->getContext().setInfoThermo(error_code); \
        } \
    } \
    catch (...) { \
        if (g_thermo) { \
            g_thermo->getContext().setInfoThermo(-999); \
        } \
    }

void TC_init() {
    try {
        if (!g_thermo) {
            g_thermo = new Thermochimica::ThermoClass();
        }
    }
    catch (...) {
        // Cannot set error code if allocation fails
        g_thermo = nullptr;
    }
}

void TC_cleanup() {
    try {
        delete g_thermo;
        g_thermo = nullptr;
    }
    catch (...) {
        // Suppress all exceptions in cleanup
        g_thermo = nullptr;
    }
}

void TC_loadDatabase(const char* filename) {
    try {
        TC_init();
        if (g_thermo && filename) {
            g_thermo->loadDatabase(filename);
        }
    }
    TC_CATCH_EXCEPTIONS(Thermochimica::ErrorCode::kDataFileNotFound)
}

// Deprecated: Use TC_loadDatabase instead
void TC_setThermoFilename(const char* filename) {
    TC_loadDatabase(filename);
}

// Deprecated: Database is loaded immediately by TC_loadDatabase
void TC_parseDataFile() {
    // No-op: loadDatabase already parses the file
    // This function is deprecated and does nothing
}

void TC_setTemperaturePressure(double T, double P) {
    try {
        TC_init();
        if (g_thermo) {
            g_thermo->setTemperaturePressure(T, P);
        }
    }
    TC_CATCH_EXCEPTIONS(Thermochimica::ErrorCode::kTemperatureOutOfRange)
}

void TC_setElementMass(int atomicNum, double mass) {
    try {
        TC_init();
        if (g_thermo) {
            g_thermo->setElementMass(atomicNum, mass);
        }
    }
    TC_CATCH_EXCEPTIONS(Thermochimica::ErrorCode::kCompositionOutOfRange)
}

void TC_solve() {
    try {
        TC_init();
        if (g_thermo) {
            g_thermo->calculate();
        }
    }
    TC_CATCH_EXCEPTIONS(Thermochimica::ErrorCode::kGEMSolverDidNotConverge)
}

int TC_getInfo() {
    try {
        if (!g_thermo) return -1;
        return g_thermo->getInfoCode();
    }
    catch (...) {
        return -999;  // Unknown exception
    }
}

double TC_getGibbsEnergy() {
    try {
        if (!g_thermo) return 0.0;
        return g_thermo->getGibbsEnergy();
    }
    catch (...) {
        return 0.0;
    }
}

void TC_reset() {
    try {
        if (g_thermo) {
            g_thermo->reset();
        }
    }
    TC_CATCH_EXCEPTIONS(-999)  // Unknown exception
}

void TC_resetAll() {
    try {
        if (g_thermo) {
            g_thermo->resetAll();
        }
    }
    TC_CATCH_EXCEPTIONS(-999)  // Unknown exception
}

} // extern "C"
