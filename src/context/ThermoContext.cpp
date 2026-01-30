#include "thermochimica/ThermoContext.hpp"

namespace Thermochimica {

ThermoContext::ThermoContext()
    : thermo(std::make_unique<ThermoState>())
    , io(std::make_unique<ThermoIO>())
    , gem(std::make_unique<GEMState>())
    , parser(std::make_unique<ParserState>())
    , submin(std::make_unique<SubMinState>())
    , ctz(std::make_unique<CTZState>())
    , reinit(std::make_unique<ReinitState>())
    , phaseConstraints(std::make_unique<PhaseConstraints>())
{
}

ThermoContext::~ThermoContext() = default;

ThermoContext::ThermoContext(ThermoContext&&) noexcept = default;

ThermoContext& ThermoContext::operator=(ThermoContext&&) noexcept = default;

void ThermoContext::resetThermo() {
    // Reset solver state but keep database loaded
    io->INFOThermo = 0;
    io->resetElementMasses();

    thermo->nConPhases = 0;
    thermo->nSolnPhases = 0;

    gem->reset();
    submin->reset();
    phaseConstraints->reset();
}

void ThermoContext::resetAll() {
    // Full reset including database
    thermo->reset();
    io->reset();
    gem->reset();
    parser->clear();
    submin->reset();
    ctz->reset();
    reinit->clear();
    phaseConstraints->clear();
}

void ThermoContext::allocateFromParser() {
    // Allocate GEM solver and reinit state arrays
    // Note: thermo->allocate() is already called in transferToThermoState
    // BEFORE data is transferred, so we don't call it again here
    int nElements = parser->nElementsCS;
    int nSolnPhases = parser->nSolnPhasesSysCS;
    int nSpecies = parser->nSpeciesCS;
    // Pure condensed phases = total species - species in solution phases
    int nConPhases = nSpecies - parser->nSpeciesPhaseCS(nSolnPhases);

    gem->allocate(nElements, nSolnPhases, nSpecies);
    reinit->allocate(nElements, nSpecies);
    phaseConstraints->allocate(nSolnPhases, nConPhases);
}

} // namespace Thermochimica
