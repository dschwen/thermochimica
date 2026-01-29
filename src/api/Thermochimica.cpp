#include "thermochimica/Thermochimica.hpp"
#include "thermochimica/parser/ChemSageParser.hpp"
#include "thermochimica/solver/GEMSolver.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>

namespace Thermochimica {

// Element symbols for lookup
static const char* kElementSymbols[] = {
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
};

void thermochimica(ThermoContext& ctx) {
    // Initialize
    init(ctx);
    if (ctx.io->INFOThermo != 0) return;

    // Check system
    checkSystem(ctx);
    if (ctx.io->INFOThermo != 0) return;

    // Compute thermodynamic data
    compThermoData(ctx);
    if (ctx.io->INFOThermo != 0) return;

    // Setup
    setup(ctx);
    if (ctx.io->INFOThermo != 0) return;

    // Solve
    solve(ctx);
}

void init(ThermoContext& ctx) {
    ctx.io->INFOThermo = 0;

    // Initialize tolerances
    ctx.thermo->tolerances.initDefaults();

    // Convert input units to internal units (K, atm, moles)
    ctx.io->dTemperature = ctx.io->convertTemperatureToKelvin(ctx.io->dTemperature);
    ctx.io->dPressure = ctx.io->convertPressureToAtm(ctx.io->dPressure);
}

void checkSystem(ThermoContext& ctx) {
    auto& io = *ctx.io;
    auto& thermo = *ctx.thermo;

    // Check temperature (must be positive, not NaN, and within bounds)
    if (std::isnan(io.dTemperature) || io.dTemperature <= 0.0 || io.dTemperature > 10000.0) {
        io.INFOThermo = ErrorCode::kTemperatureOutOfRange;
        return;
    }

    // Check pressure (must be positive, not NaN, and within bounds)
    if (std::isnan(io.dPressure) || io.dPressure <= 0.0 || io.dPressure > 1e10) {
        io.INFOThermo = ErrorCode::kPressureOutOfRange;
        return;
    }

    // Check that database is loaded
    if (thermo.nSpecies == 0) {
        io.INFOThermo = ErrorCode::kNoSpeciesInSystem;
        return;
    }

    // Check element masses - must not be NaN or negative
    for (int i = 0; i < Constants::kMaxIsotopes; ++i) {
        if (std::isnan(io.dElementMass[i]) || io.dElementMass[i] < 0.0) {
            io.INFOThermo = ErrorCode::kCompositionOutOfRange;
            return;
        }
    }

    // Check that at least one element has mass
    bool hasElement = false;
    for (int i = 0; i < Constants::kMaxIsotopes; ++i) {
        if (io.dElementMass[i] > 0.0) {
            hasElement = true;
            break;
        }
    }

    if (!hasElement) {
        io.INFOThermo = ErrorCode::kCompositionOutOfRange;
        return;
    }
}

void compThermoData(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;
    auto& parser = *ctx.parser;
    auto& io = *ctx.io;

    double T = io.dTemperature;
    double R = Constants::kIdealGasConstant;
    double logT = std::log(T);

    // Pre-compute temperature terms for standard Gibbs energy polynomial
    // G = A1 + A2*T + A3*T*ln(T) + A4*T^2 + A5*T^3 + A6/T + additional terms
    double dGibbsCoeff[6];
    dGibbsCoeff[0] = 1.0;              // Multiplies A1 (constant)
    dGibbsCoeff[1] = T;                // Multiplies A2 (linear in T)
    dGibbsCoeff[2] = T * logT;         // Multiplies A3 (T*ln(T))
    dGibbsCoeff[3] = T * T;            // Multiplies A4 (T^2)
    dGibbsCoeff[4] = T * T * T;        // Multiplies A5 (T^3)
    dGibbsCoeff[5] = 1.0 / T;          // Multiplies A6 (1/T)

    // Compute standard Gibbs energy for each species
    for (int i = 0; i < thermo.nSpecies; ++i) {
        double G = 0.0;

        // Get number of Gibbs equations for this species
        int nEqs = (i < parser.nGibbsEqSpecies.size()) ? parser.nGibbsEqSpecies(i) : 1;
        if (nEqs <= 0) nEqs = 1;

        // Find the correct temperature range equation
        int iEq = 0;
        for (int eq = 0; eq < nEqs - 1; ++eq) {
            int baseIdx = eq * ParserState::kNumGibbsCoeff;
            double Tmax = parser.dGibbsCoeffSpeciesTemp(i, baseIdx);
            if (T <= Tmax) {
                iEq = eq;
                break;
            }
            iEq = eq + 1;
        }

        // Get base index for this equation's coefficients
        int baseIdx = iEq * ParserState::kNumGibbsCoeff;

        // Compute Gibbs energy from polynomial coefficients
        // Coefficients are stored at indices 1-6 (index 0 is Tmax)
        for (int k = 0; k < 6; ++k) {
            G += parser.dGibbsCoeffSpeciesTemp(i, baseIdx + k + 1) * dGibbsCoeff[k];
        }

        // Add additional terms (up to 3 pairs of coefficient + exponent)
        // Stored at indices 7-8, 9-10, 11-12
        for (int addTerm = 0; addTerm < 3; ++addTerm) {
            int coeffIdx = baseIdx + 7 + addTerm * 2;
            int expIdx = coeffIdx + 1;

            if (expIdx < parser.dGibbsCoeffSpeciesTemp.cols()) {
                double coeff = parser.dGibbsCoeffSpeciesTemp(i, coeffIdx);
                double exp = parser.dGibbsCoeffSpeciesTemp(i, expIdx);

                if (std::abs(coeff) > 1e-30) {
                    if (std::abs(exp - 99.0) < 0.5) {
                        // Exponent 99 means ln(T)
                        G += coeff * logT;
                    } else if (std::abs(exp) > 1e-10) {
                        // General power term: coeff * T^exp
                        G += coeff * std::pow(T, exp);
                    }
                }
            }
        }

        // Store the Gibbs energy
        thermo.dStdGibbsEnergy(i) = G;

        // Normalize by RT for chemical potential
        thermo.dChemicalPotential(i) = G / (R * T);
    }

    // Compute temperature-dependent excess Gibbs energy parameters
    // G = L0 + L1*T + L2*T*ln(T) + L3*T^2 + L4*T^3 + L5/T
    for (int i = 0; i < thermo.nParam; ++i) {
        double L0 = parser.dRegularParamCS(i, 0);
        double L1 = parser.dRegularParamCS(i, 1);
        double L2 = parser.dRegularParamCS(i, 2);
        double L3 = parser.dRegularParamCS(i, 3);
        double L4 = parser.dRegularParamCS(i, 4);
        double L5 = parser.dRegularParamCS(i, 5);

        thermo.dExcessGibbsParam(i) = L0 + L1 * T + L2 * T * logT +
                                      L3 * T * T + L4 * T * T * T + L5 / T;
    }

    // Compute element amounts
    thermo.dMolesElement.setZero();
    for (int i = 0; i < Constants::kNumElementsPT; ++i) {
        if (io.dElementMass[i] > 0.0) {
            // Find element in system
            for (int j = 0; j < thermo.nElements; ++j) {
                // Match by atomic number
                int atomicNum = getAtomicNumber(thermo.cElementName[j]);
                if (atomicNum == i) {
                    thermo.dMolesElement(j) = io.dElementMass[i];
                    break;
                }
            }
        }
    }
}

void setup(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;

    // Initialize species moles
    thermo.dMolesSpecies.setZero();

    // Run leveling solver to get initial phase assemblage
    levelingSolver(ctx);

    if (ctx.io->INFOThermo != 0) return;

    // Run post-leveling solver
    postLevelingSolver(ctx);

    // Initialize GEM solver state
    GEMSolver::init(ctx);
}

void solve(ThermoContext& ctx) {
    GEMSolver::solve(ctx);

    // Post-process if successful
    if (ctx.io->INFOThermo == 0) {
        postProcess(ctx);
    }
}

void parseCSDataFile(ThermoContext& ctx, const std::string& filename) {
    setThermoFilename(ctx, filename);
    parseCSDataFile(ctx);
}

void parseCSDataFile(ThermoContext& ctx) {
    ChemSageParser::parse(ctx);
}

void setThermoFilename(ThermoContext& ctx, const std::string& filename) {
    ctx.io->cThermoFileName = filename;
}

void setTemperaturePressure(ThermoContext& ctx, double temperature, double pressure) {
    ctx.io->dTemperature = temperature;
    ctx.io->dPressure = pressure;
}

void setTemperature(ThermoContext& ctx, double temperature) {
    ctx.io->dTemperature = temperature;
}

void setPressure(ThermoContext& ctx, double pressure) {
    ctx.io->dPressure = pressure;
}

void setElementMass(ThermoContext& ctx, int atomicNumber, double mass) {
    if (atomicNumber >= 1 && atomicNumber < Constants::kMaxIsotopes) {
        ctx.io->dElementMass[atomicNumber] = mass;
    }
}

void setElementMass(ThermoContext& ctx, const std::string& elementName, double mass) {
    int atomicNum = getAtomicNumber(elementName);
    if (atomicNum > 0) {
        setElementMass(ctx, atomicNum, mass);
    }
}

void setUnits(ThermoContext& ctx,
              const std::string& tempUnit,
              const std::string& pressUnit,
              const std::string& massUnit) {
    setUnitTemperature(ctx, tempUnit);
    setUnitPressure(ctx, pressUnit);
    setUnitMass(ctx, massUnit);
}

void setUnitTemperature(ThermoContext& ctx, const std::string& unit) {
    ctx.io->cInputUnitTemperature = unit;
}

void setUnitPressure(ThermoContext& ctx, const std::string& unit) {
    ctx.io->cInputUnitPressure = unit;
}

void setUnitMass(ThermoContext& ctx, const std::string& unit) {
    ctx.io->cInputUnitMass = unit;
}

void setStandardUnits(ThermoContext& ctx) {
    setUnits(ctx, "K", "atm", "moles");
}

void setFuzzyStoich(ThermoContext& ctx, bool enable) {
    ctx.io->lFuzzyStoich = enable;
}

void setFuzzyMagnitude(ThermoContext& ctx, double magnitude) {
    ctx.io->dFuzzMag = magnitude;
}

void setPrintResultsMode(ThermoContext& ctx, int mode) {
    ctx.io->iPrintResultsMode = mode;
}

void setWriteJSON(ThermoContext& ctx, bool enable) {
    ctx.io->lWriteJSON = enable;
}

void setHeatCapacityEntropyEnthalpy(ThermoContext& ctx, bool enable) {
    ctx.io->lHeatCapacityEntropyEnthalpy = enable;
}

std::pair<double, int> getOutputChemPot(const ThermoContext& ctx,
                                        const std::string& elementName) {
    auto& thermo = *ctx.thermo;

    int elemIdx = thermo.getElementIndex(elementName);
    if (elemIdx < 0) {
        return {0.0, -1};
    }

    return {thermo.dElementPotential(elemIdx), 0};
}

std::tuple<double, double, int> getOutputSolnSpecies(
    const ThermoContext& ctx,
    const std::string& /*phaseName*/,
    const std::string& speciesName) {

    auto& thermo = *ctx.thermo;

    int speciesIdx = thermo.getSpeciesIndex(speciesName);
    if (speciesIdx < 0) {
        return {0.0, 0.0, -1};
    }

    return {thermo.dMolFraction(speciesIdx),
            thermo.dChemicalPotential(speciesIdx), 0};
}

std::tuple<double, double, int> getOutputMolSpecies(
    const ThermoContext& ctx,
    const std::string& /*phaseName*/,
    const std::string& speciesName) {

    auto& thermo = *ctx.thermo;

    int speciesIdx = thermo.getSpeciesIndex(speciesName);
    if (speciesIdx < 0) {
        return {0.0, 0.0, -1};
    }

    return {thermo.dMolesSpecies(speciesIdx),
            thermo.dMolFraction(speciesIdx), 0};
}

std::pair<double, int> getMolesPhase(const ThermoContext& ctx,
                                     const std::string& phaseName) {
    auto& thermo = *ctx.thermo;

    int phaseIdx = thermo.getPhaseIndex(phaseName);
    if (phaseIdx < 0) {
        return {0.0, -1};
    }

    // Find phase in assemblage
    for (int i = 0; i < thermo.nElements; ++i) {
        if (thermo.iAssemblage(i) == -(phaseIdx + 1)) {
            return {thermo.dMolesPhase(i), 0};
        }
    }

    return {0.0, -1};
}

std::pair<double, int> getSolnPhaseMol(const ThermoContext& ctx,
                                       const std::string& phaseName) {
    return getMolesPhase(ctx, phaseName);
}

std::pair<double, int> getPureConPhaseMol(const ThermoContext& ctx,
                                          const std::string& phaseName) {
    // For pure condensed phases, search by species name
    auto& thermo = *ctx.thermo;

    int speciesIdx = thermo.getSpeciesIndex(phaseName);
    if (speciesIdx < 0) {
        return {0.0, -1};
    }

    for (int i = 0; i < thermo.nConPhases; ++i) {
        if (thermo.iAssemblage(i) == speciesIdx + 1) {
            return {thermo.dMolesPhase(i), 0};
        }
    }

    return {0.0, -1};
}

std::pair<double, int> getElementMolesInPhase(
    const ThermoContext& ctx,
    const std::string& elementName,
    const std::string& phaseName) {

    auto& thermo = *ctx.thermo;

    int elemIdx = thermo.getElementIndex(elementName);
    int phaseIdx = thermo.getPhaseIndex(phaseName);

    if (elemIdx < 0 || phaseIdx < 0) {
        return {0.0, -1};
    }

    // Sum element contribution from all species in phase
    double sum = 0.0;
    int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
    int iLast = thermo.nSpeciesPhase(phaseIdx + 1);

    for (int i = iFirst; i < iLast; ++i) {
        sum += thermo.dMolesSpecies(i) * thermo.dStoichSpecies(i, elemIdx);
    }

    return {sum, 0};
}

int getPhaseIndex(const ThermoContext& ctx, const std::string& phaseName) {
    return ctx.thermo->getPhaseIndex(phaseName);
}

std::pair<double, int> getOutputSiteFraction(
    const ThermoContext& ctx,
    const std::string& phaseName,
    int /*sublattice*/,
    const std::string& /*constituentName*/) {

    auto& thermo = *ctx.thermo;

    int phaseIdx = thermo.getPhaseIndex(phaseName);
    if (phaseIdx < 0) {
        return {0.0, -1};
    }

    int subIdx = thermo.iPhaseSublattice(phaseIdx);
    if (subIdx <= 0 || static_cast<size_t>(subIdx - 1) >= thermo.dSiteFraction.size()) {
        return {0.0, -1};
    }

    // Find constituent index
    // Simplified - would need to search constituent names
    return {0.0, -1};
}

bool isPhaseGas(const ThermoContext& ctx, const std::string& phaseName) {
    auto& thermo = *ctx.thermo;
    int phaseIdx = thermo.getPhaseIndex(phaseName);
    if (phaseIdx < 0 || phaseIdx >= static_cast<int>(thermo.iSolnPhaseType.size())) {
        return false;
    }

    Constants::PhaseType type = thermo.iSolnPhaseType[phaseIdx];
    // IDMX (ideal mixing) is often used for gas phases
    return type == Constants::PhaseType::IDMX;
}

bool isPhaseMQM(const ThermoContext& ctx, const std::string& phaseName) {
    auto& thermo = *ctx.thermo;
    int phaseIdx = thermo.getPhaseIndex(phaseName);
    if (phaseIdx < 0 || phaseIdx >= static_cast<int>(thermo.iSolnPhaseType.size())) {
        return false;
    }

    Constants::PhaseType type = thermo.iSolnPhaseType[phaseIdx];
    return type == Constants::PhaseType::SUBG || type == Constants::PhaseType::SUBQ;
}

double getGibbsEnergy(const ThermoContext& ctx) {
    return ctx.io->dGibbsEnergySys;
}

double getHeatCapacity(const ThermoContext& ctx) {
    return ctx.io->dHeatCapacity;
}

double getEntropy(const ThermoContext& ctx) {
    return ctx.io->dEntropy;
}

double getEnthalpy(const ThermoContext& ctx) {
    return ctx.io->dEnthalpy;
}

int getNumberElementsDatabase(const ThermoContext& ctx) {
    return ctx.thermo->nElements;
}

std::string getElementAtIndex(const ThermoContext& ctx, int index) {
    if (index >= 0 && index < static_cast<int>(ctx.thermo->cElementName.size())) {
        return ctx.thermo->cElementName[index];
    }
    return "";
}

std::pair<int, int> getNumberPhasesSystem(const ThermoContext& ctx) {
    return {ctx.thermo->nSolnPhasesSys, ctx.thermo->nConPhasesSys};
}

std::string getPhaseNameAtIndex(const ThermoContext& ctx, int index) {
    if (index >= 0 && index < static_cast<int>(ctx.thermo->cSolnPhaseName.size())) {
        return ctx.thermo->cSolnPhaseName[index];
    }
    return "";
}

int getNumberSpeciesSystem(const ThermoContext& ctx) {
    return ctx.thermo->nSpecies;
}

std::string getSpeciesAtIndex(const ThermoContext& ctx, int index) {
    if (index >= 0 && index < static_cast<int>(ctx.thermo->cSpeciesName.size())) {
        return ctx.thermo->cSpeciesName[index];
    }
    return "";
}

void saveReinitData(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;
    auto& reinit = *ctx.reinit;

    reinit.saveState(thermo.iAssemblage,
                     thermo.dChemicalPotential,
                     thermo.dMolesPhase,
                     thermo.dElementPotential,
                     thermo.dMolFraction);

    ctx.io->lReinitAvailable = true;
}

void setReinitRequested(ThermoContext& ctx, bool requested) {
    ctx.io->lReinitRequested = requested;
}

bool isReinitDataAvailable(const ThermoContext& ctx) {
    return ctx.io->lReinitAvailable && ctx.reinit->isValid();
}

std::tuple<std::vector<int>, std::vector<double>, std::vector<double>,
           std::vector<double>, std::vector<double>>
getReinitData(const ThermoContext& ctx) {
    auto& reinit = *ctx.reinit;

    std::vector<int> assemblage(reinit.iAssemblage_Old.data(),
                                reinit.iAssemblage_Old.data() + reinit.iAssemblage_Old.size());
    std::vector<double> chemPot(reinit.dChemicalPotential_Old.data(),
                                reinit.dChemicalPotential_Old.data() + reinit.dChemicalPotential_Old.size());
    std::vector<double> molesPhase(reinit.dMolesPhase_Old.data(),
                                   reinit.dMolesPhase_Old.data() + reinit.dMolesPhase_Old.size());
    std::vector<double> elemPot(reinit.dElementPotential_Old.data(),
                                reinit.dElementPotential_Old.data() + reinit.dElementPotential_Old.size());
    std::vector<double> molFrac(reinit.dMolFraction_Old.data(),
                                reinit.dMolFraction_Old.data() + reinit.dMolFraction_Old.size());

    return {assemblage, chemPot, molesPhase, elemPot, molFrac};
}

void setReinitData(ThermoContext& ctx,
                   const std::vector<int>& assemblage,
                   const std::vector<double>& chemPot,
                   const std::vector<double>& molesPhase,
                   const std::vector<double>& elemPot,
                   const std::vector<double>& molFrac) {
    auto& reinit = *ctx.reinit;

    reinit.iAssemblage_Old = Eigen::Map<const Eigen::VectorXi>(assemblage.data(), assemblage.size());
    reinit.dChemicalPotential_Old = Eigen::Map<const Eigen::VectorXd>(chemPot.data(), chemPot.size());
    reinit.dMolesPhase_Old = Eigen::Map<const Eigen::VectorXd>(molesPhase.data(), molesPhase.size());
    reinit.dElementPotential_Old = Eigen::Map<const Eigen::VectorXd>(elemPot.data(), elemPot.size());
    reinit.dMolFraction_Old = Eigen::Map<const Eigen::VectorXd>(molFrac.data(), molFrac.size());

    ctx.io->lReinitAvailable = true;
}

void resetThermo(ThermoContext& ctx) {
    ctx.resetThermo();
}

void resetThermoAll(ThermoContext& ctx) {
    ctx.resetAll();
}

void printResults(ThermoContext& ctx) {
    printResultsDetailed(ctx);
}

void writeJSON(ThermoContext& ctx, bool append) {
    (void)ctx;
    (void)append;
    static bool warned = false;
    if (!warned) {
        std::cerr << "Warning: writeJSON is not yet implemented\n";
        warned = true;
    }
}

int getAtomicNumber(const std::string& symbol) {
    std::string upper = symbol;
    std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);

    // Handle special cases
    if (upper == "E-" || upper == "E(-)" || upper.substr(0, 2) == "E(") {
        return 0;  // Electron
    }

    // Search element symbols
    for (int i = 1; i <= 118; ++i) {
        std::string elemUpper = kElementSymbols[i];
        std::transform(elemUpper.begin(), elemUpper.end(), elemUpper.begin(), ::toupper);
        if (upper == elemUpper) {
            return i;
        }
    }

    return -1;  // Not found
}

std::string getElementSymbol(int atomicNumber) {
    if (atomicNumber >= 1 && atomicNumber <= 118) {
        return kElementSymbols[atomicNumber];
    }
    return "";
}

const char* getErrorMessage(int errorCode) {
    return ErrorCode::getMessage(errorCode);
}

} // namespace Thermochimica
