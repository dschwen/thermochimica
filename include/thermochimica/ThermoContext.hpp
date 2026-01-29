#pragma once

#include <memory>
#include "context/ThermoState.hpp"
#include "context/ThermoIO.hpp"
#include "context/GEMState.hpp"
#include "context/ParserState.hpp"
#include "context/SubMinState.hpp"
#include "context/CTZState.hpp"
#include "context/ReinitState.hpp"

namespace Thermochimica {

/// Main context object that replaces Fortran module globals
/// All computation state is contained within this object
/// Pass ThermoContext& to all functions for computation
class ThermoContext {
public:
    /// Core thermodynamic state (ModuleThermo.f90)
    std::unique_ptr<ThermoState> thermo;

    /// Input/Output state (ModuleThermoIO.f90)
    std::unique_ptr<ThermoIO> io;

    /// GEM solver state (ModuleGEMSolver.f90)
    std::unique_ptr<GEMState> gem;

    /// Parser temporary state (ModuleParseCS.f90)
    std::unique_ptr<ParserState> parser;

    /// Subminimization state (ModuleSubMin.f90)
    std::unique_ptr<SubMinState> submin;

    /// Common Tangent Zone state (ModuleCTZ.f90)
    std::unique_ptr<CTZState> ctz;

    /// Reinitialization state (ModuleReinit.f90)
    std::unique_ptr<ReinitState> reinit;

    /// Constructor - initializes all state objects
    ThermoContext();

    /// Destructor
    ~ThermoContext();

    /// Move constructor
    ThermoContext(ThermoContext&&) noexcept;

    /// Move assignment
    ThermoContext& operator=(ThermoContext&&) noexcept;

    // Deleted copy operations (context is not copyable)
    ThermoContext(const ThermoContext&) = delete;
    ThermoContext& operator=(const ThermoContext&) = delete;

    /// Get the current error/info code
    int infoThermo() const { return io->INFOThermo; }

    /// Set the error/info code
    void setInfoThermo(int code) { io->INFOThermo = code; }

    /// Check if computation was successful
    bool isSuccess() const { return io->INFOThermo == 0; }

    /// Reset solver state for new calculation (keeps database)
    void resetThermo();

    /// Full reset including database
    void resetAll();

    /// Check if database is loaded
    bool isDatabaseLoaded() const { return thermo->nSpecies > 0; }

    /// Allocate arrays after parsing
    void allocateFromParser();
};

} // namespace Thermochimica
