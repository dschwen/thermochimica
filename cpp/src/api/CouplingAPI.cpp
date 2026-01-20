// C-compatible API for external coupling
// Placeholder for C interface functions

#include "thermochimica/Thermochimica.hpp"

extern "C" {

// Global context for C API (singleton pattern for backward compatibility)
static Thermochimica::ThermoContext* g_ctx = nullptr;

void TC_init() {
    if (!g_ctx) {
        g_ctx = new Thermochimica::ThermoContext();
    }
}

void TC_cleanup() {
    delete g_ctx;
    g_ctx = nullptr;
}

void TC_setThermoFilename(const char* filename) {
    TC_init();
    Thermochimica::setThermoFilename(*g_ctx, filename);
}

void TC_parseDataFile() {
    TC_init();
    Thermochimica::parseCSDataFile(*g_ctx);
}

void TC_setTemperaturePressure(double T, double P) {
    TC_init();
    Thermochimica::setTemperaturePressure(*g_ctx, T, P);
}

void TC_setElementMass(int atomicNum, double mass) {
    TC_init();
    Thermochimica::setElementMass(*g_ctx, atomicNum, mass);
}

void TC_solve() {
    TC_init();
    Thermochimica::thermochimica(*g_ctx);
}

int TC_getInfo() {
    if (!g_ctx) return -1;
    return g_ctx->io->INFOThermo;
}

double TC_getGibbsEnergy() {
    if (!g_ctx) return 0.0;
    return Thermochimica::getGibbsEnergy(*g_ctx);
}

void TC_reset() {
    if (g_ctx) {
        Thermochimica::resetThermo(*g_ctx);
    }
}

void TC_resetAll() {
    if (g_ctx) {
        Thermochimica::resetThermoAll(*g_ctx);
    }
}

} // extern "C"
