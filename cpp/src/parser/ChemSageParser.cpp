#include "thermochimica/parser/ChemSageParser.hpp"
#include "thermochimica/util/ErrorCodes.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace Thermochimica {

int ChemSageParser::parse(ThermoContext& ctx, const std::string& filename) {
    ctx.io->cThermoFileName = filename;
    return parse(ctx);
}

int ChemSageParser::parse(ThermoContext& ctx) {
    ctx.io->INFOThermo = 0;
    ctx.parser->INFO = 0;
    ctx.parser->clear();

    const std::string& filename = ctx.io->cThermoFileName;

    if (filename.empty()) {
        ctx.io->INFOThermo = ErrorCode::kNoDataFileSpecified;
        return ctx.io->INFOThermo;
    }

    // Try to open the file
    std::ifstream file(filename);
    if (!file.is_open()) {
        // Try with DATA_DIRECTORY prefix
        #ifdef DATA_DIRECTORY
        std::string fullPath = std::string(DATA_DIRECTORY) + filename;
        file.open(fullPath);
        #endif

        if (!file.is_open()) {
            ctx.io->INFOThermo = ErrorCode::kDataFileNotFound;
            return ctx.io->INFOThermo;
        }
    }

    // Parse header section
    int result = parseHeader(ctx, file);
    if (result != 0) {
        ctx.io->INFOThermo = result;
        return result;
    }

    // Parse data block section
    result = parseDataBlock(ctx, file);
    if (result != 0) {
        ctx.io->INFOThermo = result;
        return result;
    }

    // Transfer parsed data to main thermo state
    transferToThermoState(ctx);

    file.close();
    return 0;
}

int ChemSageParser::parseHeader(ThermoContext& ctx, std::ifstream& file) {
    auto& p = *ctx.parser;
    std::string line;

    // Line 1: System title (ignored)
    if (!std::getline(file, line)) {
        return 101;
    }

    // Line 2: nElements, nSolnPhases, iGasPhase, [nSpeciesPerPhase...], nPureSpecies
    // Format interpretation follows Fortran ParseCSHeader.f90:
    // - First read: nElements, nSolnPhases, iGasPhase
    // - If iGasPhase == 0: no gas phase, nSolnPhases is decremented
    // - Reread line as: nElements, (dummies), nSpeciesPhaseCS(1:nSolnPhases), nSpeciesCS
    auto tokens = readTokens(file);
    if (tokens.size() < 4) {
        return 102;
    }

    p.nElementsCS = std::stoi(tokens[0]);
    p.nSolnPhasesSysCS = std::stoi(tokens[1]);
    int iGasPhase = std::stoi(tokens[2]);

    // Check limits
    if (p.nSolnPhasesSysCS > ParserState::kSolnPhasesSysMax) {
        return 8;
    }

    // Adjust for gas phase presence (follow Fortran logic)
    int nDummies = 1;  // Number of dummy values to skip when rereading
    if (iGasPhase == 0) {
        // No gas phase - number of solution phases is decremented
        p.nSolnPhasesSysCS--;
        nDummies = 1;
    } else {
        // Gas phase exists
        nDummies = 1;  // tokens[1] is a dummy
    }

    // Allocate parser arrays
    int nPhases = std::max(1, p.nSolnPhasesSysCS);
    p.nSpeciesPhaseCS.resize(nPhases + 1);
    p.nSpeciesPhaseCS.setZero();
    p.nParamPhaseCS.resize(nPhases + 1);
    p.nParamPhaseCS.setZero();
    p.nMagParamPhaseCS.resize(nPhases + 1);
    p.nMagParamPhaseCS.setZero();

    p.cElementNameCS.resize(p.nElementsCS);
    p.dAtomicMassCS.resize(p.nElementsCS);
    p.dAtomicMassCS.setZero();

    p.cSolnPhaseNameCS.resize(nPhases);
    p.cSolnPhaseTypeCS.resize(nPhases);
    p.iPhaseSublatticeCS.resize(nPhases);
    p.iPhaseSublatticeCS.setZero();

    p.nPairsSROCS.resize(nPhases, 2);
    p.nPairsSROCS.setZero();

    // Parse line 2 using Fortran-style interpretation:
    // tokens: [nElements, dummy1, nSpeciesPhase(1), ..., nSpeciesPhase(nPhases), nPureSpecies]
    // After nElements (tokens[0]), skip nDummies, then read phase species counts, then pure species
    int tokenOffset = 1 + nDummies;  // Skip nElements and dummies

    // Read species counts per phase
    if (p.nSolnPhasesSysCS > 0) {
        for (int i = 0; i < p.nSolnPhasesSysCS && tokenOffset + i < static_cast<int>(tokens.size()) - 1; ++i) {
            p.nSpeciesPhaseCS(i + 1) = std::stoi(tokens[tokenOffset + i]);
        }
    }

    // Last token is number of pure condensed species
    if (!tokens.empty()) {
        p.nSpeciesCS = std::stoi(tokens.back());
    }

    // Check for no solution phases
    if (p.nSolnPhasesSysCS == 1 && p.nSpeciesPhaseCS(1) == 0) {
        p.nSolnPhasesSysCS = 0;
    }

    // Store max species per phase
    p.nMaxSpeciesPhaseCS = p.nSpeciesPhaseCS.maxCoeff();

    // Convert to cumulative counts
    for (int i = 1; i <= p.nSolnPhasesSysCS; ++i) {
        p.nSpeciesPhaseCS(i) += p.nSpeciesPhaseCS(i - 1);
    }

    // Total species = pure + solution species
    p.nSpeciesCS += p.nSpeciesPhaseCS(p.nSolnPhasesSysCS);

    // Allocate remaining arrays
    p.dStoichSpeciesCS.resize(p.nSpeciesCS, p.nElementsCS);
    p.dStoichSpeciesCS.setZero();

    p.cSpeciesNameCS.resize(p.nSpeciesCS);
    p.nGibbsEqSpecies.resize(p.nSpeciesCS);
    p.nGibbsEqSpecies.setOnes();

    p.dGibbsCoeffSpeciesTemp.resize(p.nSpeciesCS, ParserState::kNumGibbsCoeff * ParserState::kMaxGibbsEqs);
    p.dGibbsCoeffSpeciesTemp.setZero();

    p.iPhaseCS.resize(p.nSpeciesCS);
    p.iPhaseCS.setZero();

    p.iParticlesPerMoleCS.resize(p.nSpeciesCS);
    p.iParticlesPerMoleCS.setOnes();

    p.dQKTOParamsCS.resize(p.nSpeciesCS, 2);
    p.dQKTOParamsCS.setZero();

    // Line 3: Element names
    tokens = readTokens(file);
    if (static_cast<int>(tokens.size()) < p.nElementsCS) {
        return 103;
    }

    for (int i = 0; i < p.nElementsCS; ++i) {
        p.cElementNameCS[i] = tokens[i];
        // Convert electron notation
        if (p.cElementNameCS[i].substr(0, 2) == "e(") {
            p.cElementNameCS[i] = "e-";
        }
    }

    // Allocate sublattice arrays
    p.nSublatticePhaseCS.resize(p.nSolnPhasesSysCS);
    p.nSublatticePhaseCS.setZero();

    p.dStoichSublatticeCS.resize(p.nSolnPhasesSysCS, ParserState::kMaxSublattice);
    p.dStoichSublatticeCS.setZero();

    p.nConstituentSublatticeCS.resize(p.nSolnPhasesSysCS, ParserState::kMaxSublattice);
    p.nConstituentSublatticeCS.setZero();

    p.dZetaSpeciesCS.resize(p.nSolnPhasesSysCS, p.nMaxSpeciesPhaseCS);
    p.dZetaSpeciesCS.setZero();

    p.nInterpolationOverrideCS.resize(p.nSolnPhasesSysCS);
    p.nInterpolationOverrideCS.setZero();

    // Line 4: Atomic masses
    auto masses = readDoubles(file, p.nElementsCS);
    if (static_cast<int>(masses.size()) < p.nElementsCS) {
        return 104;
    }

    for (int i = 0; i < p.nElementsCS; ++i) {
        p.dAtomicMassCS(i) = masses[i];
    }

    // Line 5: Temperature dependence terms (2 lines, verification only)
    auto tempTerms1 = readInts(file, 7);
    auto tempTerms2 = readInts(file, 7);

    // Verify standard temperature dependence
    std::vector<int> expected = {6, 1, 2, 3, 4, 5, 6};
    bool valid = true;
    for (int i = 0; i < 7 && i < static_cast<int>(tempTerms2.size()); ++i) {
        if (tempTerms2[i] != expected[i]) {
            valid = false;
            break;
        }
    }

    if (!valid) {
        return 105;
    }

    return 0;
}

int ChemSageParser::parseDataBlock(ThermoContext& ctx, std::ifstream& file) {
    auto& p = *ctx.parser;
    std::string line;

    int speciesIndex = 0;

    // Parse solution phases
    for (int iPhase = 0; iPhase < p.nSolnPhasesSysCS; ++iPhase) {
        // Read phase name (Entry 1 - separate line, Fortran uses FORMAT(A25))
        if (!std::getline(file, line)) {
            return 1100 + iPhase;
        }
        // Trim whitespace
        size_t start = line.find_first_not_of(" \t");
        if (start == std::string::npos) start = 0;
        size_t end = line.find_last_not_of(" \t");
        p.cSolnPhaseNameCS[iPhase] = line.substr(start, end - start + 1);

        // Read phase type (Entry 2 - separate line)
        auto tokens = readTokens(file);
        if (tokens.empty()) {
            return 1200 + iPhase;
        }
        p.cSolnPhaseTypeCS[iPhase] = tokens[0];

        // Check if phase type is supported
        if (!ParserState::isPhaseTypeSupported(p.cSolnPhaseTypeCS[iPhase])) {
            return ErrorCode::kUnsupportedPhaseType;
        }

        // Get phase type
        Constants::PhaseType phaseType = ParserState::getPhaseType(p.cSolnPhaseTypeCS[iPhase]);

        // Parse phase-specific data based on type
        int nSpeciesInPhase = p.nSpeciesPhaseCS(iPhase + 1) - p.nSpeciesPhaseCS(iPhase);

        // For sublattice phases, read additional structure data
        if (phaseType == Constants::PhaseType::SUBL ||
            phaseType == Constants::PhaseType::SUBLM ||
            phaseType == Constants::PhaseType::SUBI ||
            phaseType == Constants::PhaseType::SUBM) {
            int result = parseSUBLPhase(ctx, file, iPhase);
            if (result != 0) return result;
        } else if (phaseType == Constants::PhaseType::SUBG ||
                   phaseType == Constants::PhaseType::SUBQ) {
            int result = parseSUBGPhase(ctx, file, iPhase);
            if (result != 0) return result;
        }

        // Parse each species in this phase
        for (int i = 0; i < nSpeciesInPhase; ++i) {
            int result = parseSolutionPhase(ctx, file, speciesIndex);
            if (result != 0) return result;

            p.iPhaseCS(speciesIndex) = iPhase + 1;  // 1-based phase index
            ++speciesIndex;
        }

        // Parse mixing parameters for this phase (only for non-ideal phases)
        // IDMX (ideal mixing) has no mixing parameters
        if (phaseType != Constants::PhaseType::IDMX) {
            tokens = readTokens(file);
            if (!tokens.empty()) {
                try {
                    int nParams = std::stoi(tokens[0]);
                    if (nParams > 0) {
                        int result = parseMixingParameters(ctx, file, iPhase, nParams);
                        if (result != 0) return result;
                    }
                } catch (const std::exception& e) {
                    return 1700 + iPhase;  // Mixing params conversion error
                }
            }
        }
    }

    // Parse pure condensed phases
    int nPureSpecies = p.nSpeciesCS - p.nSpeciesPhaseCS(p.nSolnPhasesSysCS);
    for (int i = 0; i < nPureSpecies; ++i) {
        int result = parseSolutionPhase(ctx, file, speciesIndex);
        if (result != 0) return result;

        p.iPhaseCS(speciesIndex) = 0;  // Pure condensed phase
        ++speciesIndex;
    }

    return 0;
}

int ChemSageParser::parseSolutionPhase(ThermoContext& ctx, std::ifstream& file,
                                       int speciesIndex) {
    auto& p = *ctx.parser;
    std::string line;

    // Entry 3: Read species name (Fortran FORMAT A26)
    if (!std::getline(file, line)) {
        return 1300 + speciesIndex;
    }
    // Trim whitespace
    size_t start = line.find_first_not_of(" \t");
    if (start == std::string::npos) {
        p.cSpeciesNameCS[speciesIndex] = "";
    } else {
        size_t end = line.find_last_not_of(" \t");
        p.cSpeciesNameCS[speciesIndex] = line.substr(start, end - start + 1);
    }

    // Check for particles per mole (e.g., "):2" means 2 particles per mole)
    size_t colonPos = p.cSpeciesNameCS[speciesIndex].rfind(':');
    if (colonPos != std::string::npos && colonPos >= 2) {
        std::string prefix = p.cSpeciesNameCS[speciesIndex].substr(colonPos - 2, 2);
        if (prefix == "):") {
            try {
                p.iParticlesPerMoleCS(speciesIndex) =
                    std::stoi(p.cSpeciesNameCS[speciesIndex].substr(colonPos + 1));
            } catch (...) {
                p.iParticlesPerMoleCS(speciesIndex) = 1;
            }
        }
    }

    // Entry 4: Read iGibbsEqType, nGibbsEqSpecies, stoichiometry (all on one line)
    auto tokens = readTokens(file);
    if (static_cast<int>(tokens.size()) < 2 + p.nElementsCS) {
        return 1400 + speciesIndex;
    }

    int iGibbsEqType = 0;
    int nGibbsEqs = 0;
    try {
        iGibbsEqType = std::stoi(tokens[0]);
        nGibbsEqs = std::stoi(tokens[1]);
    } catch (const std::exception& e) {
        return 1401 + speciesIndex;  // Conversion error
    }
    p.nGibbsEqSpecies(speciesIndex) = nGibbsEqs;

    // Check supported Gibbs equation types (4, 16, 1, 13)
    if (!(iGibbsEqType == 4 || iGibbsEqType == 16 || iGibbsEqType == 1 || iGibbsEqType == 13)) {
        return 1402 + speciesIndex;  // Unsupported equation type
    }

    // Read stoichiometry
    try {
        for (int j = 0; j < p.nElementsCS && 2 + j < static_cast<int>(tokens.size()); ++j) {
            p.dStoichSpeciesCS(speciesIndex, j) = std::stod(tokens[2 + j]) *
                static_cast<double>(p.iParticlesPerMoleCS(speciesIndex));
        }
    } catch (const std::exception& e) {
        return 1403 + speciesIndex;  // Stoichiometry conversion error
    }

    // Parse Gibbs energy equations
    int result = parseGibbsEquations(ctx, file, speciesIndex, nGibbsEqs, iGibbsEqType);
    if (result != 0) return result;

    return 0;
}

int ChemSageParser::parseGibbsEquations(ThermoContext& ctx, std::ifstream& file,
                                        int speciesIndex, int nGibbsEqs, int iGibbsEqType) {
    auto& p = *ctx.parser;

    for (int iEq = 0; iEq < nGibbsEqs; ++iEq) {
        // Entry 5: Read first 7 values (Tmax + 6 standard Gibbs energy coefficients)
        auto values = readDoubles(file, 7);
        if (values.size() < 7) {
            return 1500 + speciesIndex;
        }

        int baseIdx = iEq * ParserState::kNumGibbsCoeff;

        // Store Tmax and 6 coefficients (indices 0-6)
        for (int k = 0; k < 7 && baseIdx + k < p.dGibbsCoeffSpeciesTemp.cols(); ++k) {
            p.dGibbsCoeffSpeciesTemp(speciesIndex, baseIdx + k) = values[k];
        }

        // For types 4 and 16, there's an additional line with more terms
        if (iGibbsEqType == 4 || iGibbsEqType == 16) {
            // Read the number of additional terms (first value of next line)
            auto addTerms = readTokens(file);
            if (addTerms.empty()) {
                continue;  // No additional terms
            }

            try {
                int nAddTerms = std::stoi(addTerms[0]);

                // Read additional coefficients based on nAddTerms
                if (nAddTerms >= 1 && addTerms.size() >= 3) {
                    p.dGibbsCoeffSpeciesTemp(speciesIndex, baseIdx + 7) = std::stod(addTerms[1]);
                    p.dGibbsCoeffSpeciesTemp(speciesIndex, baseIdx + 8) = std::stod(addTerms[2]);
                }
                if (nAddTerms >= 2 && addTerms.size() >= 5) {
                    p.dGibbsCoeffSpeciesTemp(speciesIndex, baseIdx + 9) = std::stod(addTerms[3]);
                    p.dGibbsCoeffSpeciesTemp(speciesIndex, baseIdx + 10) = std::stod(addTerms[4]);
                }
                if (nAddTerms >= 3 && addTerms.size() >= 7) {
                    p.dGibbsCoeffSpeciesTemp(speciesIndex, baseIdx + 11) = std::stod(addTerms[5]);
                    p.dGibbsCoeffSpeciesTemp(speciesIndex, baseIdx + 12) = std::stod(addTerms[6]);
                }
            } catch (const std::exception& e) {
                return 1600 + speciesIndex;  // Additional terms conversion error
            }
        }
    }

    return 0;
}

int ChemSageParser::parseSUBLPhase(ThermoContext& ctx, std::ifstream& file,
                                   int phaseIndex) {
    auto& p = *ctx.parser;

    // Read number of sublattices
    auto tokens = readTokens(file);
    if (tokens.empty()) {
        return 1500 + phaseIndex;
    }

    int nSublattices = std::stoi(tokens[0]);
    p.nSublatticePhaseCS(phaseIndex) = nSublattices;
    p.iPhaseSublatticeCS(phaseIndex) = ++p.nCountSublatticeCS;

    // Read sublattice stoichiometry
    auto stoich = readDoubles(file, nSublattices);
    for (int i = 0; i < nSublattices && i < static_cast<int>(stoich.size()); ++i) {
        p.dStoichSublatticeCS(phaseIndex, i) = stoich[i];
    }

    // Read constituents per sublattice
    tokens = readTokens(file);
    for (int i = 0; i < nSublattices && i < static_cast<int>(tokens.size()); ++i) {
        p.nConstituentSublatticeCS(phaseIndex, i) = std::stoi(tokens[i]);
    }

    return 0;
}

int ChemSageParser::parseSUBGPhase(ThermoContext& ctx, std::ifstream& file,
                                   int phaseIndex) {
    auto& p = *ctx.parser;

    // Read quadruplet information
    auto tokens = readTokens(file);
    if (tokens.size() < 2) {
        return 1600 + phaseIndex;
    }

    p.nPairsSROCS(phaseIndex, 0) = std::stoi(tokens[0]);
    p.nPairsSROCS(phaseIndex, 1) = std::stoi(tokens[1]);
    p.iPhaseSublatticeCS(phaseIndex) = ++p.nCountSublatticeCS;

    return 0;
}

int ChemSageParser::parseMixingParameters(ThermoContext& ctx, std::ifstream& file,
                                          int phaseIndex, int nParams) {
    auto& p = *ctx.parser;

    for (int i = 0; i < nParams; ++i) {
        // Read parameter line
        auto tokens = readTokens(file);
        // Store parameter data
        // Implementation depends on specific phase type
        ++p.nParamCS;
    }

    p.nParamPhaseCS(phaseIndex + 1) = p.nParamCS;

    return 0;
}

void ChemSageParser::transferToThermoState(ThermoContext& ctx) {
    auto& p = *ctx.parser;
    auto& t = *ctx.thermo;

    // Allocate main thermo state
    t.allocate(p.nElementsCS, p.nSpeciesCS, p.nSolnPhasesSysCS, p.nParamCS);

    // Transfer dimensions
    t.nElements = p.nElementsCS;
    t.nSpecies = p.nSpeciesCS;
    t.nSolnPhasesSys = p.nSolnPhasesSysCS;
    t.nConPhasesSys = p.nSpeciesCS - p.nSpeciesPhaseCS(p.nSolnPhasesSysCS);
    t.nParam = p.nParamCS;
    t.nCountSublattice = p.nCountSublatticeCS;

    // Transfer element names and atomic masses
    for (int i = 0; i < p.nElementsCS; ++i) {
        t.cElementName[i] = p.cElementNameCS[i];
        t.dAtomicMass(i) = p.dAtomicMassCS(i);
    }

    // Transfer species data
    for (int i = 0; i < p.nSpeciesCS; ++i) {
        t.cSpeciesName[i] = p.cSpeciesNameCS[i];
        t.iPhase(i) = p.iPhaseCS(i);
        t.iParticlesPerMole(i) = p.iParticlesPerMoleCS(i);

        for (int j = 0; j < p.nElementsCS; ++j) {
            t.dStoichSpecies(i, j) = p.dStoichSpeciesCS(i, j);
        }
    }

    // Transfer phase data
    for (int i = 0; i < p.nSolnPhasesSysCS; ++i) {
        t.cSolnPhaseName[i] = p.cSolnPhaseNameCS[i];
        t.cSolnPhaseType[i] = p.cSolnPhaseTypeCS[i];
        t.nSpeciesPhase(i + 1) = p.nSpeciesPhaseCS(i + 1);
        t.iPhaseSublattice(i) = p.iPhaseSublatticeCS(i);
    }

    // Allocate context arrays
    ctx.allocateFromParser();
}

std::vector<std::string> ChemSageParser::readTokens(std::ifstream& file) {
    std::vector<std::string> tokens;
    std::string line;

    if (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        while (iss >> token) {
            tokens.push_back(token);
        }
    }

    return tokens;
}

std::vector<double> ChemSageParser::readDoubles(std::ifstream& file, int count) {
    std::vector<double> values;
    std::string line;

    while (static_cast<int>(values.size()) < count && std::getline(file, line)) {
        std::istringstream iss(line);
        double val;
        while (iss >> val && static_cast<int>(values.size()) < count) {
            values.push_back(val);
        }
    }

    return values;
}

std::vector<int> ChemSageParser::readInts(std::ifstream& file, int count) {
    std::vector<int> values;
    std::string line;

    while (static_cast<int>(values.size()) < count && std::getline(file, line)) {
        std::istringstream iss(line);
        int val;
        while (iss >> val && static_cast<int>(values.size()) < count) {
            values.push_back(val);
        }
    }

    return values;
}

std::string ChemSageParser::trim(const std::string& str) {
    auto start = str.find_first_not_of(" \t\n\r");
    if (start == std::string::npos) return "";
    auto end = str.find_last_not_of(" \t\n\r");
    return str.substr(start, end - start + 1);
}

} // namespace Thermochimica
