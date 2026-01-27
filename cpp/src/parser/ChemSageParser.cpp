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
        file = std::ifstream(fullPath);  // Create fresh ifstream
        #endif

        if (!file.is_open()) {
            ctx.io->INFOThermo = ErrorCode::kDataFileNotFound;
            return ctx.io->INFOThermo;
        }
    }

    try {
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
    } catch (const std::exception& e) {
        // Catch any uncaught parsing exceptions
        ctx.io->INFOThermo = ErrorCode::kParserError;
        return ErrorCode::kParserError;
    }

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

    // Line 2 (possibly spanning multiple lines): nElements, nSolnPhases, iGasPhase, [nSpeciesPerPhase...], nPureSpecies
    // Format interpretation follows Fortran ParseCSHeader.f90:
    // - First read: nElements, nSolnPhases, iGasPhase
    // - If iGasPhase == 0: no gas phase, nSolnPhases is decremented
    // - Reread line as: nElements, (dummies), nSpeciesPhaseCS(1:nSolnPhases), nSpeciesCS

    // Read header integers (may span multiple lines)
    // Keep reading lines of integers until we get a line that starts with a non-integer
    std::vector<std::string> tokens;
    std::streampos lastPos = file.tellg();
    while (file) {
        auto newTokens = readTokens(file);

        // Check if this line contains numeric tokens
        bool allNumeric = !newTokens.empty();
        for (const auto& tok : newTokens) {
            try {
                std::stoi(tok);
            } catch (...) {
                allNumeric = false;
                break;
            }
        }

        if (!allNumeric) {
            // This line has non-numeric tokens (element names) - rewind and stop
            file.seekg(lastPos);
            break;
        }

        tokens.insert(tokens.end(), newTokens.begin(), newTokens.end());
        lastPos = file.tellg();

        // Safety limit to prevent infinite loop
        if (tokens.size() > 200) break;
    }
    if (tokens.size() < 4) {
        return 102;
    }

    try {
        p.nElementsCS = std::stoi(tokens[0]);
        p.nSolnPhasesInFile = std::stoi(tokens[1]);  // Store original count
        p.nSolnPhasesSysCS = p.nSolnPhasesInFile;
        p.iGasPhase = std::stoi(tokens[2]);
    } catch (const std::exception& e) {
        // Header parsing failed - token is not a valid integer
        return 103;
    }

    // Check limits
    if (p.nSolnPhasesSysCS > ParserState::kSolnPhasesSysMax) {
        return 8;
    }

    // Adjust for gas phase presence (follow Fortran logic exactly)
    // The iGasPhase value determines how many dummy values to skip:
    // - iGasPhase == 0: no gas phase, skip 2 dummies (tokens[1] and tokens[2])
    // - iGasPhase != 0: gas phase exists, skip 1 dummy (tokens[1])
    int nDummies;
    if (p.iGasPhase == 0) {
        // No gas phase - number of solution phases is decremented
        p.nSolnPhasesSysCS--;
        nDummies = 2;  // Skip tokens[1] (nSolnPhases) and tokens[2] (iGasPhase=0)
    } else {
        // Gas phase exists
        nDummies = 1;  // Skip tokens[1] (nSolnPhases, which is a dummy when gas exists)
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
    // tokens: [nElements, nSolnPhases, iGasPhase, nSpeciesPhase(1), ..., nSpeciesPhase(nPhases), nPureSpecies]
    // After nElements (tokens[0]), skip nDummies to get to species counts
    // nDummies=1 when gas exists: skip tokens[1]
    // nDummies=2 when no gas: skip tokens[1] and tokens[2]
    int tokenOffset = 1 + nDummies;  // Skip nElements and dummies

    // Read species counts per phase
    try {
        if (p.nSolnPhasesSysCS > 0) {
            for (int i = 0; i < p.nSolnPhasesSysCS && tokenOffset + i < static_cast<int>(tokens.size()) - 1; ++i) {
                p.nSpeciesPhaseCS(i + 1) = std::stoi(tokens[tokenOffset + i]);
            }
        }

        // Last token is number of pure condensed species
        if (!tokens.empty()) {
            p.nSpeciesCS = std::stoi(tokens.back());
        }
    } catch (const std::exception& e) {
        // Species count parsing failed
        return 104;
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

    // Allocate mixing parameter arrays (estimate max params = 10 per phase)
    int maxParams = p.nSolnPhasesSysCS * 20;  // Conservative estimate
    p.iRegularParamCS.resize(maxParams, 12);  // Up to 12 columns for indices/exponents
    p.iRegularParamCS.setZero();
    p.dRegularParamCS.resize(maxParams, 6);   // 6 L coefficients
    p.dRegularParamCS.setZero();

    // Line 3: Element names - may span multiple lines
    tokens = readTokens(file, p.nElementsCS);
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
    int realPhaseIndex = 0;  // Index for actual stored phases

    // Parse solution phases
    // When iGasPhase == 0, the file contains nSolnPhasesSysCS phases (not nSolnPhasesInFile)
    // because the "absent gas phase" is not in the data, just the header count was adjusted
    int nPhasesInFile = p.nSolnPhasesSysCS;
    for (int iFilePhase = 0; iFilePhase < nPhasesInFile; ++iFilePhase) {
        // Read phase name (Entry 1 - separate line, Fortran uses FORMAT(A25))
        std::streampos posBeforeName = file.tellg();
        if (!std::getline(file, line)) {
            return 1100 + iFilePhase;
        }
        // Trim whitespace
        size_t start = line.find_first_not_of(" \t");
        std::string phaseName;
        if (start == std::string::npos) {
            phaseName = "";
        } else {
            size_t end = line.find_last_not_of(" \t");
            phaseName = line.substr(start, end - start + 1);
        }

        // Read phase type (Entry 2 - separate line)
        auto tokens = readTokens(file);
        if (tokens.empty()) {
            return 1200 + iFilePhase;
        }
        std::string phaseType = tokens[0];

        // Check if this is an empty placeholder phase
        // Fortran skips: 1) phases where first 6 chars are spaces, 2) phase name exactly "Delete"
        // For " Delete" (with leading space), it's NOT skipped - it's a real phase
        bool isEmptyPhase = (line.length() >= 6 && line.substr(0, 6) == "      ");
        // Also check for exact "Delete" at start of line (no leading whitespace)
        if (!isEmptyPhase && line.length() >= 6 && line.substr(0, 6) == "Delete") {
            isEmptyPhase = true;
        }

        if (isEmptyPhase) {
            // Empty placeholder phase - just name and type, no species
            continue;
        }

        // This is a real phase to process
        p.cSolnPhaseNameCS[realPhaseIndex] = phaseName;
        p.cSolnPhaseTypeCS[realPhaseIndex] = phaseType;


        // Check if phase type is supported
        if (!ParserState::isPhaseTypeSupported(phaseType)) {
            return ErrorCode::kUnsupportedPhaseType;
        }

        // Get phase type enum
        Constants::PhaseType phaseTypeEnum = ParserState::getPhaseType(phaseType);

        // Parse phase-specific data based on type
        int nSpeciesInPhase = p.nSpeciesPhaseCS(realPhaseIndex + 1) - p.nSpeciesPhaseCS(realPhaseIndex);

        // For sublattice phases, read additional structure data
        // Note: SUBI has a different format - species are given directly without sublattice header
        // Some SUBL files also have post-species structure like SUBI
        bool sublHasPreHeader = false;
        if (phaseTypeEnum == Constants::PhaseType::SUBL ||
            phaseTypeEnum == Constants::PhaseType::SUBLM ||
            phaseTypeEnum == Constants::PhaseType::SUBM) {
            int result = parseSUBLPhase(ctx, file, realPhaseIndex, sublHasPreHeader);
            if (result != 0) return result;
        } else if (phaseTypeEnum == Constants::PhaseType::SUBI) {
            // SUBI phases don't have the same sublattice header - mark as having sublattice data
            // but don't try to parse extra header lines
            p.iPhaseSublatticeCS(realPhaseIndex) = ++p.nCountSublatticeCS;
        } else if (phaseTypeEnum == Constants::PhaseType::SUBG ||
                   phaseTypeEnum == Constants::PhaseType::SUBQ) {
            int result = parseSUBGPhase(ctx, file, realPhaseIndex);
            if (result != 0) return result;
        } else if (phaseTypeEnum == Constants::PhaseType::RKMPM) {
            // RKMPM phases have a sublattice stoichiometry header line before species
            // Format: stoich1 stoich2 (e.g., "0.333333        0.280000")
            auto stoichValues = readDoubles(file, 2);
            if (stoichValues.size() >= 2) {
                p.nSublatticePhaseCS(realPhaseIndex) = 2;  // RKMPM has 2 sublattices
                p.dStoichSublatticeCS(realPhaseIndex, 0) = stoichValues[0];
                p.dStoichSublatticeCS(realPhaseIndex, 1) = stoichValues[1];
            }
        }

        // Parse each species in this phase
        for (int i = 0; i < nSpeciesInPhase; ++i) {
            int result = parseSolutionPhase(ctx, file, speciesIndex, phaseTypeEnum);
            if (result != 0) {
                return result;
            }

            p.iPhaseCS(speciesIndex) = realPhaseIndex + 1;  // 1-based phase index
            ++speciesIndex;
        }

        // For SUBI phases and SUBL/SUBLM/SUBM phases (regardless of pre-header),
        // parse the sublattice structure data that comes AFTER species
        bool needsPostSpeciesData = (phaseTypeEnum == Constants::PhaseType::SUBI) ||
            (phaseTypeEnum == Constants::PhaseType::SUBL) ||
            (phaseTypeEnum == Constants::PhaseType::SUBLM) ||
            (phaseTypeEnum == Constants::PhaseType::SUBM);

        if (needsPostSpeciesData) {
            int result = parseSUBIExcessData(ctx, file, realPhaseIndex);
            if (result != 0) return result;
        }
        // Parse mixing parameters for this phase (only for non-ideal phases)
        // IDMX (ideal mixing) has no mixing parameters
        else if (phaseTypeEnum != Constants::PhaseType::IDMX) {
            // QKTO/RKMP use a different format: read until "0" terminator
            // Each parameter starts with type indicator (2=binary, 3=ternary)
            int result = parseMixingParametersLoop(ctx, file, realPhaseIndex, phaseTypeEnum);
            if (result != 0) return result;
        }

        ++realPhaseIndex;
    }

    // Parse pure condensed phases (no QKTO parameters)
    int nPureSpecies = p.nSpeciesCS - p.nSpeciesPhaseCS(p.nSolnPhasesSysCS);
    for (int i = 0; i < nPureSpecies; ++i) {
        int result = parseSolutionPhase(ctx, file, speciesIndex, Constants::PhaseType::Unknown);
        if (result != 0) return result;

        p.iPhaseCS(speciesIndex) = 0;  // Pure condensed phase
        ++speciesIndex;
    }

    return 0;
}

int ChemSageParser::parseSolutionPhase(ThermoContext& ctx, std::ifstream& file,
                                       int speciesIndex, Constants::PhaseType phaseType) {
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
    if (result != 0) {
        return result;
    }

    // For QKTO phases, read K and p parameters after Gibbs equations
    // These appear as a line like "1.00000  1" after each species' Gibbs data
    if (phaseType == Constants::PhaseType::QKTO) {
        auto qktoParams = readDoubles(file, 2);
        if (qktoParams.size() >= 2) {
            p.dQKTOParamsCS(speciesIndex, 0) = qktoParams[0];  // K parameter
            p.dQKTOParamsCS(speciesIndex, 1) = qktoParams[1];  // p parameter
        }
    }

    // For SUBQ/SUBG phases, read coordination numbers and charge after Gibbs equations
    // Line 1: 5 coordination numbers (or nElements + 1 values)
    // Line 2: Charge/electronegativity
    if (phaseType == Constants::PhaseType::SUBQ || phaseType == Constants::PhaseType::SUBG) {
        // Read coordination numbers (typically 5 values: z_A, z_B, z_C, z_D, z_E)
        auto coordNums = readDoubles(file, p.nElementsCS + 1);
        // Store coordination numbers if needed (in p.dCoordNumbersCS or similar)

        // Read charge/electronegativity
        auto charge = readDoubles(file, 1);
        // Store charge if needed (in p.dSpeciesChargeCS or similar)
    }

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

        // For types 4 and 16, there's an additional line after EACH temp range
        // with "nAddTerms coef1 exp1 coef2 exp2..."
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

    // For type 13, read ONE continuation line with 2 values AFTER all temp ranges
    if (iGibbsEqType == 13) {
        auto addTerms = readDoubles(file, 2);
        // Store in first temp range's extra slots if needed
        if (addTerms.size() >= 2) {
            p.dGibbsCoeffSpeciesTemp(speciesIndex, 7) = addTerms[0];
            p.dGibbsCoeffSpeciesTemp(speciesIndex, 8) = addTerms[1];
        }
    }

    // For type 16 (magnetic), read magnetic ordering parameters after all Gibbs equations
    // Format: Curie/Neel temperature, magnetic moment
    if (iGibbsEqType == 16) {
        auto magData = readDoubles(file, 2);
        // Store if needed: magData[0] = Tc (Curie/Neel), magData[1] = magnetic moment
        if (magData.size() >= 2) {
            // Could store in p.dMagneticMoment and p.dCurieTemp if arrays exist
        }
    }

    return 0;
}

int ChemSageParser::parseSUBLPhase(ThermoContext& ctx, std::ifstream& file,
                                   int phaseIndex, bool& hasPreHeader) {
    auto& p = *ctx.parser;
    hasPreHeader = false;

    try {
        // Save file position to potentially seek back
        std::streampos startPos = file.tellg();

        // Read next line to check format
        std::string line;
        if (!std::getline(file, line)) {
            return 1500 + phaseIndex;
        }

        // Determine if this is a sublattice stoichiometry line (all numbers)
        // or a species name line (contains letters/colon)
        std::istringstream iss(line);
        std::vector<double> stoichValues;
        double val;
        bool isNumericLine = true;

        // Try parsing all tokens as doubles
        std::string token;
        while (iss >> token) {
            try {
                size_t pos;
                val = std::stod(token, &pos);
                if (pos != token.length()) {
                    // Not fully parsed as number
                    isNumericLine = false;
                    break;
                }
                stoichValues.push_back(val);
            } catch (...) {
                isNumericLine = false;
                break;
            }
        }

        if (!isNumericLine || stoichValues.empty()) {
            // No sublattice stoichiometry header - this is a species name
            // Seek back so species parser can read this line
            file.seekg(startPos);

            // Mark that there's no pre-header (post-species structure will be needed)
            hasPreHeader = false;
            p.nSublatticePhaseCS(phaseIndex) = 2;  // Default assumption
            p.iPhaseSublatticeCS(phaseIndex) = ++p.nCountSublatticeCS;
            return 0;
        }

        // Has stoichiometry header - number of values = number of sublattices
        hasPreHeader = true;
        int nSublattices = static_cast<int>(stoichValues.size());
        p.nSublatticePhaseCS(phaseIndex) = nSublattices;
        p.iPhaseSublatticeCS(phaseIndex) = ++p.nCountSublatticeCS;

        // Store the stoichiometry values we already parsed
        for (int i = 0; i < nSublattices; ++i) {
            p.dStoichSublatticeCS(phaseIndex, i) = stoichValues[i];
        }

        // Note: In this format, constituents per sublattice is determined from species
        // The format just has stoichiometry, then species directly
        // No separate "constituents per sublattice" line before species
    } catch (const std::exception& e) {
        return 1501 + phaseIndex;  // SUBL parsing error
    }

    return 0;
}

int ChemSageParser::parseSUBGPhase(ThermoContext& ctx, std::ifstream& file,
                                   int phaseIndex) {
    auto& p = *ctx.parser;

    try {
        // Read quadruplet information
        auto tokens = readTokens(file);
        if (tokens.size() < 2) {
            return 1600 + phaseIndex;
        }

        p.nPairsSROCS(phaseIndex, 0) = std::stoi(tokens[0]);
        p.nPairsSROCS(phaseIndex, 1) = std::stoi(tokens[1]);
        p.iPhaseSublatticeCS(phaseIndex) = ++p.nCountSublatticeCS;
    } catch (const std::exception& e) {
        return 1601 + phaseIndex;  // SUBG parsing error
    }

    return 0;
}

int ChemSageParser::parseSUBIExcessData(ThermoContext& ctx, std::ifstream& file,
                                         int phaseIndex) {
    auto& p = *ctx.parser;
    std::string line;

    // SUBI phases have sublattice structure and excess parameters AFTER the species
    // Format (based on CsTe-1.dat):
    //   nSublattice nConstituent
    //   constituent names (one per line for first sublattice)
    //   constituent names (one line for second sublattice)
    //   number of species on sublattice 2
    //   charge values
    //   stoichiometry
    //   iConstituentSublattice indices
    //   more indices
    //   nParams
    //   param indices
    //   param values (6 per line)
    //   0 (terminator)

    try {
        // Read first line to determine format
        auto tokens = readTokens(file);

        if (tokens.empty()) {
            return 1700 + phaseIndex;
        }

        int nSublattice = std::stoi(tokens[0]);

        // Determine format: SUBI has (nSublattice, nConstituent), SUBL has nSublattice alone
        bool isSUBIFormat = (tokens.size() >= 2);

        if (isSUBIFormat) {
            // SUBI format: nSublattice nConstituent on same line
            int nConstituent = std::stoi(tokens[1]);

            p.nSublatticePhaseCS(phaseIndex) = nSublattice;

            // Read constituent names - one line per sublattice
            for (int i = 0; i < nSublattice; ++i) {
                if (!std::getline(file, line)) {
                    return 1701 + phaseIndex;
                }
            }

            // Read additional SUBI-specific structure data
            tokens = readTokens(file);  // nSpecies
            tokens = readTokens(file);  // charge
            tokens = readTokens(file);  // stoichiometry
            tokens = readTokens(file);  // indices line 1
            tokens = readTokens(file);  // indices line 2
        } else {
            // SUBL format: nSublattice alone, then more structure lines

            p.nSublatticePhaseCS(phaseIndex) = nSublattice;

            // Read stoichiometry (nSublattice values)
            tokens = readTokens(file);

            // Read constituents per sublattice (nSublattice values)
            tokens = readTokens(file);

            // Read constituent names - one line per sublattice
            for (int i = 0; i < nSublattice; ++i) {
                if (!std::getline(file, line)) {
                    return 1701 + phaseIndex;
                }
            }

            // Read iConstituentSublattice indices (2 lines)
            tokens = readTokens(file);
            tokens = readTokens(file);
        }

        // Now read excess parameters in a loop
        // For SUBLM: TWO sections - magnetic params (ending in 0), then regular excess (ending in 0)
        // For SUBI/SUBL: ONE section of excess params (ending in 0)
        // Format: type (2=binary, 3=ternary), indices line, value lines (count from last index)

        // Read first terminator section (magnetic params for SUBLM, excess for SUBI/SUBL)
        int blockCount = 0;
        int sectionCount = 0;
        while (true) {
            tokens = readTokens(file);
            if (tokens.empty()) {
                return 1704 + phaseIndex;
            }

            int paramType = std::stoi(tokens[0]);
            if (paramType == 0) {
                // Found terminator for this section
                ++sectionCount;
                blockCount = 0;

                // For SUBLM, there's a second section (regular excess) after magnetic
                // Check if next line starts another param section or is another terminator
                std::streampos pos = file.tellg();
                auto peekTokens = readTokens(file);
                if (!peekTokens.empty()) {
                    try {
                        int nextVal = std::stoi(peekTokens[0]);
                        // If it's a small int (2 or 3), it's a param type - continue parsing
                        if (nextVal >= 2 && nextVal <= 3) {
                            file.seekg(pos);  // Seek back to re-read this line
                            continue;  // Parse second section
                        }
                        // If it's another 0, this is the second terminator - consume it and exit
                        if (nextVal == 0) {
                            // Don't seek back - we consumed the second terminator
                            break;
                        }
                    } catch (...) {
                        // Not an int - must be next phase name
                    }
                }
                // Seek back and exit - done with all sections
                file.seekg(pos);
                break;
            }

            // Read indices line
            tokens = readTokens(file);

            // The last index indicates the number of value lines (order/multiplicity)
            int nValueLines = 1;
            if (!tokens.empty()) {
                try {
                    nValueLines = std::stoi(tokens.back());
                    if (nValueLines < 1) nValueLines = 1;
                    if (nValueLines > 10) nValueLines = 1;  // Sanity check
                } catch (...) {
                    nValueLines = 1;
                }
            }

            // Read value lines (6 coefficients each)
            for (int i = 0; i < nValueLines; ++i) {
                auto values = readDoubles(file, 6);
            }
            ++blockCount;
        }

    } catch (const std::exception& e) {
        // If parsing fails, try to recover by skipping to terminator
        while (std::getline(file, line)) {
            size_t start = line.find_first_not_of(" \t");
            if (start != std::string::npos && line[start] == '0' &&
                (start + 1 >= line.length() || line.find_first_not_of(" \t", start + 1) == std::string::npos)) {
                // Found standalone "0" - terminator
                break;
            }
        }
    }

    return 0;
}

int ChemSageParser::parseMixingParametersLoop(ThermoContext& ctx, std::ifstream& file,
                                               int phaseIndex, Constants::PhaseType phaseType) {
    auto& p = *ctx.parser;

    // Mixing parameters format varies by phase type
    // For QKTO: each parameter has type (2=binary, 3=ternary), then indices+values
    // For RKMP: nParams count, then each param has 3 indices + 6 values
    // Continue reading until type=0 is encountered

    while (true) {
        // Read first token - either parameter type (QKTO) or count (RKMP)
        auto tokens = readTokens(file);
        if (tokens.empty()) {
            break;
        }

        int firstVal = 0;
        try {
            firstVal = std::stoi(tokens[0]);
        } catch (...) {
            break;
        }

        // Check for terminator (0)
        if (firstVal == 0) {
            break;  // End of mixing parameters
        }

        // For QKTO: firstVal is the parameter type (2 or 3)
        // For type 2 (binary): 4 indices (i, j, p, q) + 6 L coefficients
        // For type 3 (ternary): 6 indices (i, j, k, p, q, r) + 6 L coefficients
        if (phaseType == Constants::PhaseType::QKTO) {
            int nSpeciesParam = firstVal;  // 2 for binary, 3 for ternary
            int nIndices = (nSpeciesParam == 2) ? 4 : 6;

            // Read indices + 6 L coefficient values
            auto data = readDoubles(file, nIndices + 6);
            if (static_cast<int>(data.size()) < nIndices + 6) {
                break;  // Incomplete data
            }

            // Store in iRegularParamCS if we have space
            int iParam = p.nParamCS;
            if (iParam < p.iRegularParamCS.rows()) {
                // Column 1: number of species in parameter (2 or 3)
                p.iRegularParamCS(iParam, 1) = nSpeciesParam;

                // Columns 2,3,(4): species indices (convert 1-based file index to 0-based)
                p.iRegularParamCS(iParam, 2) = static_cast<int>(data[0]) - 1;  // species a (0-based)
                p.iRegularParamCS(iParam, 3) = static_cast<int>(data[1]) - 1;  // species b (0-based)
                if (nSpeciesParam == 3) {
                    p.iRegularParamCS(iParam, 4) = static_cast<int>(data[2]) - 1;  // species c (0-based)
                }

                // Exponents p, q, (r) at columns nSpeciesParam+2, +3, (+4)
                int expOffset = (nSpeciesParam == 2) ? 2 : 3;  // After species indices
                p.iRegularParamCS(iParam, nSpeciesParam + 2) = static_cast<int>(data[expOffset]);      // p
                p.iRegularParamCS(iParam, nSpeciesParam + 3) = static_cast<int>(data[expOffset + 1]);  // q
                if (nSpeciesParam == 3) {
                    p.iRegularParamCS(iParam, nSpeciesParam + 4) = static_cast<int>(data[expOffset + 2]);  // r
                }

                // Store L coefficients (6 values for temperature polynomial)
                p.dRegularParamCS(iParam, 0) = data[nIndices];      // L0
                p.dRegularParamCS(iParam, 1) = data[nIndices + 1];  // L1
                p.dRegularParamCS(iParam, 2) = data[nIndices + 2];  // L2
                p.dRegularParamCS(iParam, 3) = data[nIndices + 3];  // L3
                p.dRegularParamCS(iParam, 4) = data[nIndices + 4];  // L4
                p.dRegularParamCS(iParam, 5) = data[nIndices + 5];  // L5
            }

            ++p.nParamCS;
        }
        // RKMPM format: similar to QKTO but with separate lines
        // Type (2 or 3) on one line, indices on next, then value lines
        else if (phaseType == Constants::PhaseType::RKMPM) {
            int paramType = firstVal;  // 2=binary, 3=ternary

            // Read indices (i, j, k where k is order/exponent)
            auto indices = readTokens(file);
            if (indices.size() < 3) break;

            int lastIndex = std::stoi(indices.back());  // Order/exponent
            int nValueLines = std::max(1, lastIndex);  // Number of value lines

            // Read value lines
            for (int v = 0; v < nValueLines; ++v) {
                auto values = readDoubles(file, 6);
                // Store if needed
            }

            ++p.nParamCS;
            // Continue loop to read more params until 0
        }
        // For RKMP and other types: firstVal is nParams count
        else {
            // RKMP format: firstVal is parameter count
            for (int i = 0; i < firstVal; ++i) {
                // Read indices (typically 3)
                auto indices = readTokens(file);
                if (indices.size() < 3) break;

                int iParam = p.nParamCS;
                if (iParam < p.iRegularParamCS.rows()) {
                    // Store indices
                    p.iRegularParamCS(iParam, 1) = 2;  // Binary
                    p.iRegularParamCS(iParam, 2) = std::stoi(indices[0]);
                    p.iRegularParamCS(iParam, 3) = std::stoi(indices[1]);
                    p.iRegularParamCS(iParam, 4) = std::stoi(indices[2]);
                }

                // Read 6 coefficient values
                auto values = readDoubles(file, 6);
                if (values.size() < 6) break;

                if (iParam < p.dRegularParamCS.rows()) {
                    for (int j = 0; j < 6; ++j) {
                        p.dRegularParamCS(iParam, j) = values[j];
                    }
                }

                ++p.nParamCS;
            }
            // After RKMP params, look for magnetic params (0 or count)
            auto magTokens = readTokens(file);
            if (!magTokens.empty()) {
                try {
                    int nMagParams = std::stoi(magTokens[0]);
                    for (int i = 0; i < nMagParams; ++i) {
                        readTokens(file);  // indices
                        readDoubles(file, 6);  // values
                    }
                } catch (...) {}
            }
            break;  // RKMP has explicit count, so we're done
        }
    }

    p.nParamPhaseCS(phaseIndex + 1) = p.nParamCS;

    return 0;
}

int ChemSageParser::parseMixingParameters(ThermoContext& ctx, std::ifstream& file,
                                          int phaseIndex, int nParams) {
    auto& p = *ctx.parser;


    // For RKMP, each parameter has: 3 indices (i, j, k) + 6 values = 9 total values
    // The indices are on one line, values on the next (possibly spanning lines)
    int paramsRead = 0;
    for (int i = 0; i < nParams; ++i) {
        // Read indices line (3 integers)
        auto indices = readTokens(file);

        // Check if this looks like mixing param indices (should be 3 integers)
        // or if it's actually the magnetic params line (1 small integer)
        if (indices.size() == 1) {
            // This might be magnetic params count, not mixing indices
            // Check if it's a small integer (magnetic params are typically 0 or small)
            try {
                int val = std::stoi(indices[0]);
                if (val <= 10) {  // Likely magnetic params count, not mixing index
                    // Handle magnetic params
                    for (int j = 0; j < val; ++j) {
                        readTokens(file);  // indices
                        readDoubles(file, 6);  // values
                    }
                    break;  // Exit mixing params loop
                }
            } catch (...) {}
        }

        if (indices.size() < 3) {
            break;
        }

        // Read 6 coefficient values
        auto values = readDoubles(file, 6);

        if (values.size() < 6) {
            break;
        }

        ++paramsRead;
        ++p.nParamCS;
    }

    // If we read all expected params, look for magnetic params line
    if (paramsRead == nParams) {
        auto magTokens = readTokens(file);
        if (!magTokens.empty()) {
            try {
                int nMagParams = std::stoi(magTokens[0]);
                for (int i = 0; i < nMagParams; ++i) {
                    readTokens(file);  // indices
                    readDoubles(file, 6);  // values
                }
            } catch (...) {}
        }
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

        // Transfer QKTO parameters (K, p, group)
        t.dQKTOParams(i, 0) = p.dQKTOParamsCS(i, 0);  // K parameter
        t.dQKTOParams(i, 1) = p.dQKTOParamsCS(i, 1);  // p parameter
        t.dQKTOParams(i, 2) = 0;  // Group index (default 0 - all species in same group)
    }

    // Transfer phase data
    for (int i = 0; i < p.nSolnPhasesSysCS; ++i) {
        t.cSolnPhaseName[i] = p.cSolnPhaseNameCS[i];
        t.cSolnPhaseType[i] = p.cSolnPhaseTypeCS[i];
        t.nSpeciesPhase(i + 1) = p.nSpeciesPhaseCS(i + 1);
        t.iPhaseSublattice(i) = p.iPhaseSublatticeCS(i);
        t.nParamPhase(i + 1) = p.nParamPhaseCS(i + 1);  // Cumulative param count
    }

    // Transfer mixing parameters
    for (int i = 0; i < p.nParamCS; ++i) {
        // Transfer integer param data (indices, exponents)
        for (int j = 0; j < std::min(static_cast<int>(p.iRegularParamCS.cols()),
                                     static_cast<int>(t.iRegularParam.cols())); ++j) {
            t.iRegularParam(i, j) = p.iRegularParamCS(i, j);
        }
        // Transfer L coefficients - compute temperature-dependent value at T=298.15K as initial
        // G = L0 + L1*T + L2*T*ln(T) + L3*T^2 + L4*T^3 + L5/T
        // For now, just store L0 as the base value
        t.dExcessGibbsParam(i) = p.dRegularParamCS(i, 0);
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

std::vector<std::string> ChemSageParser::readTokens(std::ifstream& file, int count) {
    std::vector<std::string> tokens;
    std::string line;

    while (static_cast<int>(tokens.size()) < count && std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        while (iss >> token && static_cast<int>(tokens.size()) < count) {
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
