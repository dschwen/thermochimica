// C-compatible API for external coupling
// Provides C interface functions for use from C code or other languages via FFI

#include "thermochimica/ThermoClass.hpp"

extern "C" {

// Global ThermoClass instance for C API (singleton pattern for backward compatibility)
static Thermochimica::ThermoClass* g_thermo = nullptr;

void TC_init() {
    if (!g_thermo) {
        g_thermo = new Thermochimica::ThermoClass();
    }
}

void TC_cleanup() {
    delete g_thermo;
    g_thermo = nullptr;
}

void TC_loadDatabase(const char* filename) {
    TC_init();
    g_thermo->loadDatabase(filename);
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
    TC_init();
    g_thermo->setTemperaturePressure(T, P);
}

void TC_setElementMass(int atomicNum, double mass) {
    TC_init();
    g_thermo->setElementMass(atomicNum, mass);
}

void TC_solve() {
    TC_init();
    g_thermo->calculate();
}

int TC_getInfo() {
    if (!g_thermo) return -1;
    return g_thermo->getInfoCode();
}

double TC_getGibbsEnergy() {
    if (!g_thermo) return 0.0;
    return g_thermo->getGibbsEnergy();
}

void TC_reset() {
    if (g_thermo) {
        g_thermo->reset();
    }
}

void TC_resetAll() {
    if (g_thermo) {
        g_thermo->resetAll();
    }
}

} // extern "C"
