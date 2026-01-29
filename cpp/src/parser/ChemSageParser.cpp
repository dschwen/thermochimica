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

    // Adjust for gas phase presence (follow Fortran ParseCSHeader.f90 logic exactly)
    // The Fortran does:
    //   1. Read nElements, nSolnPhases, iGasPhase
    //   2. If iGasPhase == 0: no gas, decrement nSolnPhases, set iGasPhase = 1
    //   3. Else: gas exists, set iGasPhase = 0
    //   4. Re-read: nElements, iDummy(1:iGasPhase+1), nSpeciesPhase(1:nSolnPhases), nSpecies
    //
    // So iDummy count is (adjusted iGasPhase + 1):
    // - No gas (original iGasPhase==0 -> adjusted=1): iDummy(1:2) = 2 dummies
    // - Gas exists (original iGasPhase!=0 -> adjusted=0): iDummy(1:1) = 1 dummy
    int adjustedGasPhase;
    if (p.iGasPhase == 0) {
        // No gas phase - number of solution phases is decremented
        p.nSolnPhasesSysCS--;
        adjustedGasPhase = 1;  // After Fortran adjustment
    } else {
        // Gas phase exists
        adjustedGasPhase = 0;  // After Fortran adjustment
    }
    int nDummies = adjustedGasPhase + 1;

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
    // Fortran re-reads as: nElements, iDummy(1:adjustedGasPhase+1), nSpeciesPhase(1:nSolnPhases), nSpecies
    // When gas exists (adjustedGasPhase=0): skip 1 dummy (tokens[1])
    // When no gas (adjustedGasPhase=1): skip 2 dummies (tokens[1], tokens[2])
    // So tokenOffset = 1 + nDummies where nDummies = adjustedGasPhase + 1
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

    // Allocate magnetic parameter arrays (for RKMPM phases)
    p.iMagneticParamCS.resize(maxParams, 7);  // Up to 7 columns: type + species indices + exponent
    p.iMagneticParamCS.setZero();
    p.dMagneticParamCS.resize(maxParams, 2);  // 2 coefficients per magnetic param
    p.dMagneticParamCS.setZero();

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

    p.nInterpolationOverrideCS.resize(p.nSolnPhasesSysCS + 1);  // +1 for 1-based indexing
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

        // Determine actual number of species to parse
        // For SUBQ/SUBG, only nPairsSROCS[0] species are actual species data
        // The rest of nSpeciesInPhase are pair fraction slots (not parsed as species)
        int nSpeciesToParse = nSpeciesInPhase;
        if (phaseTypeEnum == Constants::PhaseType::SUBQ ||
            phaseTypeEnum == Constants::PhaseType::SUBG) {
            nSpeciesToParse = p.nPairsSROCS(realPhaseIndex, 0);
        }

        // Parse each species in this phase
        for (int i = 0; i < nSpeciesToParse; ++i) {
            int result = parseSolutionPhase(ctx, file, speciesIndex, phaseTypeEnum,
                                            realPhaseIndex, i);
            if (result != 0) {
                return result;
            }

            p.iPhaseCS(speciesIndex) = realPhaseIndex + 1;  // 1-based phase index
            ++speciesIndex;
        }

        // For SUBQ/SUBG phases, need to advance speciesIndex for remaining pair slots
        // (they're allocated but not read as species)
        if (phaseTypeEnum == Constants::PhaseType::SUBQ ||
            phaseTypeEnum == Constants::PhaseType::SUBG) {
            int remainingSlots = nSpeciesInPhase - nSpeciesToParse;
            for (int i = 0; i < remainingSlots; ++i) {
                p.iPhaseCS(speciesIndex) = realPhaseIndex + 1;
                ++speciesIndex;
            }
        }

        // For SUBQ/SUBG phases, parse sublattice structure and mixing parameters
        if (phaseTypeEnum == Constants::PhaseType::SUBQ ||
            phaseTypeEnum == Constants::PhaseType::SUBG) {
            int result = parseSUBGExcessData(ctx, file, realPhaseIndex, phaseTypeEnum);
            if (result != 0) return result;
        }
        // For SUBI phases, SUBL/SUBLM/SUBM phases:
        // Parse the sublattice structure data that comes AFTER species
        // Note: SUBL/SUBLM/SUBM may have a pre-header stoichiometry line BUT
        // they ALSO have post-species structure (constituent names, indices, etc.)
        // before the mixing parameters
        else if ((phaseTypeEnum == Constants::PhaseType::SUBI) ||
                 (phaseTypeEnum == Constants::PhaseType::SUBL) ||
                 (phaseTypeEnum == Constants::PhaseType::SUBLM) ||
                 (phaseTypeEnum == Constants::PhaseType::SUBM)) {
            int result = parseSUBIExcessData(ctx, file, realPhaseIndex, phaseTypeEnum);
            if (result != 0) return result;
        }
        // Parse mixing parameters for this phase (only for non-ideal phases)
        // IDMX (ideal mixing) has no mixing parameters
        else if (phaseTypeEnum != Constants::PhaseType::IDMX) {
            // RKMPM phases have magnetic mixing parameters BEFORE regular mixing parameters
            // Both sections are terminated by "0"
            if (phaseTypeEnum == Constants::PhaseType::RKMPM) {
                int result = parseMagneticParameters(ctx, file, realPhaseIndex);
                if (result != 0) return result;
            }
            // QKTO/RKMP/RKMPM use a different format: read until "0" terminator
            // Each parameter starts with type indicator (2=binary, 3=ternary)
            int result = parseMixingParametersLoop(ctx, file, realPhaseIndex, phaseTypeEnum);
            if (result != 0) return result;
        }

        ++realPhaseIndex;
    }

    // Parse pure condensed phases (no QKTO parameters)
    int nPureSpecies = p.nSpeciesCS - p.nSpeciesPhaseCS(p.nSolnPhasesSysCS);
    for (int i = 0; i < nPureSpecies; ++i) {
        int result = parseSolutionPhase(ctx, file, speciesIndex, Constants::PhaseType::Unknown,
                                        -1, i);  // -1 = no phase, i = local index
        if (result != 0) return result;

        p.iPhaseCS(speciesIndex) = 0;  // Pure condensed phase
        ++speciesIndex;
    }

    return 0;
}

int ChemSageParser::parseSolutionPhase(ThermoContext& ctx, std::ifstream& file,
                                       int speciesIndex, Constants::PhaseType phaseType,
                                       int phaseIndex, int localSpeciesIndex) {
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

    // For SUBQ/SUBG phases, read constituent coefficients (and zeta for SUBQ) after Gibbs equations
    // Line 1: 5 constituent coefficients (z_A, z_B, z_C, z_D, z_E)
    // Line 2 (SUBQ only): per-species zeta value
    if (phaseType == Constants::PhaseType::SUBQ || phaseType == Constants::PhaseType::SUBG) {
        // Read 5 constituent coefficients (always 5, regardless of nElements)
        auto coeffs = readDoubles(file, 5);
        if (coeffs.size() < 5) {
            return 1500 + speciesIndex;
        }

        // Get the sublattice phase index (use the phase's sublattice index)
        int sublatticePhaseIdx = p.iPhaseSublatticeCS(phaseIndex) - 1;  // Convert to 0-based
        if (sublatticePhaseIdx >= 0) {
            // Ensure the 3D array is sized properly
            if (static_cast<int>(p.dConstituentCoefficientsCS.size()) <= sublatticePhaseIdx) {
                p.dConstituentCoefficientsCS.resize(sublatticePhaseIdx + 1);
            }
            if (p.dConstituentCoefficientsCS[sublatticePhaseIdx].rows() == 0) {
                p.dConstituentCoefficientsCS[sublatticePhaseIdx].resize(p.nMaxSpeciesPhaseCS, 5);
                p.dConstituentCoefficientsCS[sublatticePhaseIdx].setZero();
            }

            // Store constituent coefficients using the passed local species index
            for (int c = 0; c < 5; ++c) {
                p.dConstituentCoefficientsCS[sublatticePhaseIdx](localSpeciesIndex, c) = coeffs[c];
            }
        }

        // For SUBQ only: read per-species zeta value
        if (phaseType == Constants::PhaseType::SUBQ) {
            auto zeta = readDoubles(file, 1);
            if (zeta.empty()) {
                return 1501 + speciesIndex;
            }

            // Store zeta in dZetaSpeciesCS (indexed by [sublattice_phase, local_species])
            if (sublatticePhaseIdx >= 0 && sublatticePhaseIdx < p.dZetaSpeciesCS.rows()) {
                if (localSpeciesIndex < p.dZetaSpeciesCS.cols()) {
                    p.dZetaSpeciesCS(sublatticePhaseIdx, localSpeciesIndex) = zeta[0];
                }
            }
        }
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

    // For types 16 and 13, read magnetic/continuation line AFTER all Gibbs equations
    // For solution phases: 2 values (dGibbsMagneticCS(1:2))
    // For pure condensed phases: 4 values (dGibbsMagneticCS(1:4))
    // Note: We can't easily tell pure vs solution here, so we read 2 values
    if (iGibbsEqType == 16 || iGibbsEqType == 13) {
        auto magData = readDoubles(file, 2);
        // Store magnetic contribution: Curie/Neel temp and magnetic moment/beta
        // Could store in p.dGibbsMagneticCS if the array exists
        (void)magData;  // Suppress unused variable warning
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

int ChemSageParser::parseSUBGExcessData(ThermoContext& ctx, std::ifstream& file,
                                         int phaseIndex, Constants::PhaseType phaseType) {
    auto& p = *ctx.parser;
    std::string line;

    // SUBG/SUBQ phases have sublattice structure and mixing parameters AFTER species
    // Format (based on ParseCSDataBlockSUBG.f90):
    //   nConstituentSublatticeCS(1:2) - # constituents per sublattice
    //   constituent names (3 per line, 25 chars each) for each sublattice
    //   charges for sublattice 1
    //   chemical groups for sublattice 1
    //   charges for sublattice 2
    //   chemical groups for sublattice 2
    //   iConstituentSublattice IDs for sublattice 1
    //   iConstituentSublattice IDs for sublattice 2
    //   coordination number data
    //   mixing parameters (with Q/G prefix lines)

    try {
        int sublatticeIdx = p.iPhaseSublatticeCS(phaseIndex) - 1;  // 0-based
        if (sublatticeIdx < 0) sublatticeIdx = 0;

        // Line 1: Number of constituents per sublattice (2 values)
        auto tokens = readTokens(file);
        if (tokens.size() < 2) {
            return 1602 + phaseIndex;
        }

        int nConst1 = std::stoi(tokens[0]);
        int nConst2 = std::stoi(tokens[1]);

        // Ensure arrays are allocated
        if (static_cast<int>(p.nConstituentSublatticeCS.rows()) <= sublatticeIdx) {
            // Resize needed but difficult - just use what we have
        }
        p.nConstituentSublatticeCS(sublatticeIdx, 0) = nConst1;
        p.nConstituentSublatticeCS(sublatticeIdx, 1) = nConst2;
        p.nSublatticePhaseCS(phaseIndex) = 2;  // SUBG/SUBQ always have 2 sublattices

        // Read constituent names for sublattice 1 (3 names per line, 25 chars each)
        int nLines1 = (nConst1 + 2) / 3;  // Ceiling division
        for (int i = 0; i < nLines1; ++i) {
            if (!std::getline(file, line)) return 1603 + phaseIndex;
            // Parse 3 names of 25 chars each (store if needed)
        }

        // Read constituent names for sublattice 2
        int nLines2 = (nConst2 + 2) / 3;
        for (int i = 0; i < nLines2; ++i) {
            if (!std::getline(file, line)) return 1604 + phaseIndex;
        }

        // Read charges for sublattice 1
        readDoubles(file, nConst1);

        // Read chemical groups for sublattice 1
        readTokens(file);

        // Read charges for sublattice 2
        readDoubles(file, nConst2);

        // Read chemical groups for sublattice 2
        readTokens(file);

        // Number of pairs for iConstituentSublattice
        int nPairs = p.nPairsSROCS(phaseIndex, 0);  // Use actual species count

        // Read iConstituentSublattice IDs for sublattice 1 (nPairs values)
        readTokens(file);

        // Read iConstituentSublattice IDs for sublattice 2 (nPairs values)
        readTokens(file);

        // Read coordination number data
        // The number of coordination lines is nPairsSROCS(phaseIndex, 1)
        // Each line: j k x y za zb zx zy (4 ints + 4 doubles)
        int nCoordLines = p.nPairsSROCS(phaseIndex, 1);
        for (int i = 0; i < nCoordLines; ++i) {
            // Each line has 4 integers and 4 doubles (8 values total)
            if (!std::getline(file, line)) break;
            // Just consume the line - we'd need to store coordination data properly
        }

        // Parse mixing parameters (based on ParseCSDataBlockSUBG.f90 lines 373-416)
        // Format:
        //   Line 1: parameter type (3=ternary, 4=quaternary, or 0/negative to end)
        //   Line 2: Q/G/R/B char + 8 indices
        //   Lines 3-4: 2 ignored chi lines (6 doubles each)
        //   Line 5: 2 integers + 6 coefficients

        while (true) {
            tokens = readTokens(file);
            if (tokens.empty()) break;

            int paramType = std::stoi(tokens[0]);

            // Check for terminator (0 or negative)
            // Negative value indicates number of interpolation override lines
            if (paramType <= 0) {
                // Skip interpolation override lines if any
                int nOverrides = -paramType;
                for (int ovr = 0; ovr < nOverrides; ++ovr) {
                    if (!std::getline(file, line)) break;
                }
                break;
            }

            // Validate parameter type (3=ternary or 4=quaternary)
            if (paramType != 3 && paramType != 4) {
                // Unknown parameter type - error
                return 1610 + phaseIndex;
            }

            int iParam = p.nParamCS;

            // Read Q/G/R/B line with 8 indices
            // Format: cRegularParamCS(nParamCS), iRegularParamCS(nParamCS,2:9)
            tokens = readTokens(file);
            if (tokens.empty()) break;

            // First token is Q/G/R/B character
            char paramChar = tokens[0][0];
            if (paramChar != 'Q' && paramChar != 'G' && paramChar != 'R' && paramChar != 'B') {
                return 1611 + phaseIndex;
            }

            // Store parameter type and character
            if (iParam < p.iRegularParamCS.rows()) {
                p.iRegularParamCS(iParam, 0) = paramType;  // Type in column 0
                p.iRegularParamCS(iParam, 1) = paramType;  // Also in column 1 for consistency
                p.cRegularParamCS.push_back(paramChar);

                // Store 8 indices from tokens[1..8] into iRegularParamCS columns 2-9
                for (int idx = 0; idx < 8 && idx + 1 < static_cast<int>(tokens.size()); ++idx) {
                    p.iRegularParamCS(iParam, 2 + idx) = std::stoi(tokens[idx + 1]);
                }
            }

            // Read 2 ignored chi interpolation lines (6 doubles each)
            readDoubles(file, 6);
            readDoubles(file, 6);

            // Read coefficient line: 2 integers + 6 doubles (may span 2 lines)
            // Format: iRegularParamCS(nParamCS,10:11), dRegularParamCS(nParamCS,1:6)
            tokens = readTokens(file, 8);
            if (tokens.size() < 8) {
                // Need at least 2 ints + 6 doubles = 8 values
                return 1612 + phaseIndex;
            }

            if (iParam < p.iRegularParamCS.rows()) {
                // Store the 2 integers in columns 10-11
                p.iRegularParamCS(iParam, 10) = std::stoi(tokens[0]);
                p.iRegularParamCS(iParam, 11) = std::stoi(tokens[1]);

                // Store 6 coefficients
                for (int c = 0; c < 6; ++c) {
                    p.dRegularParamCS(iParam, c) = std::stod(tokens[2 + c]);
                }
            }

            ++p.nParamCS;
        }

        p.nParamPhaseCS(phaseIndex + 1) = p.nParamCS;

    } catch (const std::exception& e) {
        return 1605 + phaseIndex;
    }

    return 0;
}

int ChemSageParser::parseSUBIExcessData(ThermoContext& ctx, std::ifstream& file,
                                         int phaseIndex, Constants::PhaseType phaseType) {
    auto& p = *ctx.parser;
    std::string line;
    bool isSUBLM = (phaseType == Constants::PhaseType::SUBLM);

    // SUBI/SUBL phases have sublattice structure and excess parameters AFTER the species
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
            // SUBL format (ParseCSDataBlockSUBL.f90):
            // 1. nSublattice (already read)
            // 2. stoichiometry (nSublattice values)
            // 3. constituents per sublattice (nSublattice values)
            // 4. For each sublattice: ceil(nConstituents/3) lines of names (3 per line, 25 chars)
            // 5. For each sublattice: 1 line of iConstituentSublattice indices (nSpecies values)

            p.nSublatticePhaseCS(phaseIndex) = nSublattice;

            // Read stoichiometry (nSublattice values)
            tokens = readTokens(file);

            // Read constituents per sublattice (nSublattice values)
            auto constTokens = readTokens(file);
            std::vector<int> nConstituents(nSublattice, 0);
            for (int i = 0; i < nSublattice && i < static_cast<int>(constTokens.size()); ++i) {
                nConstituents[i] = std::stoi(constTokens[i]);
            }

            // Read constituent names - ceil(nConstituents[s]/3) lines per sublattice
            for (int s = 0; s < nSublattice; ++s) {
                int nLines = (nConstituents[s] + 2) / 3;  // Ceiling division
                for (int ln = 0; ln < nLines; ++ln) {
                    if (!std::getline(file, line)) {
                        return 1701 + phaseIndex;
                    }
                }
            }

            // Read iConstituentSublattice indices - one line per sublattice
            for (int s = 0; s < nSublattice; ++s) {
                tokens = readTokens(file);
            }
        }

        // Now read excess parameters in a loop
        // For SUBLM: TWO sections - magnetic params (ending in 0), then regular excess (ending in 0)
        // For SUBI/SUBL: ONE section of excess params (ending in 0)
        // Format: type (2=binary, 3=ternary), indices line, value lines (count from last index)

        // Read mixing params sections
        // For SUBLM: first section is magnetic (2 coefficients), second is regular (6 coefficients)
        // For SUBI/SUBL: only regular params (6 coefficients)
        int sectionCount = 0;
        bool inMagneticSection = isSUBLM;  // SUBLM starts with magnetic params

        while (true) {
            tokens = readTokens(file);
            if (tokens.empty()) {
                return 1704 + phaseIndex;
            }

            int paramType = std::stoi(tokens[0]);
            if (paramType == 0) {
                // Found terminator for this section
                ++sectionCount;

                if (isSUBLM && sectionCount == 1) {
                    // Just finished magnetic section, now read regular params
                    inMagneticSection = false;
                    continue;  // Parse second section (regular)
                }
                // Either non-SUBLM phase or finished both sections for SUBLM
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

            // Read value lines
            // Magnetic params: 2 coefficients per line
            // Regular params: 6 coefficients per line
            int nCoeffs = inMagneticSection ? 2 : 6;
            for (int i = 0; i < nValueLines; ++i) {
                auto values = readDoubles(file, nCoeffs);
            }
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
        // RKMP/RKMPM format: paramType (2/3/4) on one line, indices on next, then k value lines
        // Loop continues until "0" terminator is read
        // For binary (2): indices are (species1, species2, k), then k lines of 6 coefficients
        // For ternary (3): indices are (species1, species2, species3, k), then k lines
        // For quaternary (4): indices are (species1, species2, species3, species4, k), then k lines
        else if (phaseType == Constants::PhaseType::RKMP ||
                 phaseType == Constants::PhaseType::RKMPM) {
            int paramType = firstVal;  // 2=binary, 3=ternary, 4=quaternary

            // Determine number of indices to read based on param type
            int nIndices = paramType + 1;  // binary=3, ternary=4, quaternary=5

            // Read indices line
            auto indices = readTokens(file);
            if (static_cast<int>(indices.size()) < nIndices) break;

            // Parse species indices (1-based in file)
            std::vector<int> speciesIdx(paramType);
            for (int s = 0; s < paramType; ++s) {
                speciesIdx[s] = std::stoi(indices[s]);
            }

            // Last index is k (number of coefficient sets / polynomial order)
            int k = std::stoi(indices[paramType]);
            k = std::max(1, k);

            // Read k sets of 6 coefficients, each becomes a separate parameter
            for (int v = 0; v < k; ++v) {
                auto values = readDoubles(file, 6);
                if (values.size() < 6) break;

                int iParam = p.nParamCS;
                if (iParam < p.iRegularParamCS.rows()) {
                    // Column 1: number of species in parameter
                    p.iRegularParamCS(iParam, 1) = paramType;

                    // Store species indices (keep 1-based as Fortran does)
                    for (int s = 0; s < paramType; ++s) {
                        p.iRegularParamCS(iParam, 2 + s) = speciesIdx[s];
                    }

                    // Set the exponent (0, 1, 2, ... for each coefficient set)
                    p.iRegularParamCS(iParam, 2 + paramType) = v;

                    // Store L coefficients
                    for (int c = 0; c < 6; ++c) {
                        p.dRegularParamCS(iParam, c) = values[c];
                    }
                }

                ++p.nParamCS;
            }
            // Continue loop to read more params until 0
        }
        // For other phase types, skip unknown mixing format
        else {
            // Unknown phase type - try to skip to terminator
            break;
        }
    }

    p.nParamPhaseCS(phaseIndex + 1) = p.nParamCS;

    return 0;
}

int ChemSageParser::parseMagneticParameters(ThermoContext& ctx, std::ifstream& file,
                                            int phaseIndex) {
    auto& p = *ctx.parser;

    // RKMPM phases have magnetic mixing parameters that come BEFORE regular mixing parameters
    // Format is similar to regular mixing params but with only 2 coefficients per line
    // Loop until "0" terminator is read

    while (true) {
        // Read first token - parameter type (2=binary, 3=ternary, 4=quaternary) or 0 to terminate
        auto tokens = readTokens(file);
        if (tokens.empty()) {
            break;
        }

        int paramType = 0;
        try {
            paramType = std::stoi(tokens[0]);
        } catch (...) {
            break;
        }

        // Check for terminator (0)
        if (paramType == 0) {
            break;  // End of magnetic parameters
        }

        // Determine number of indices to read based on param type
        // binary (2): species1, species2, k -> 3 indices
        // ternary (3): species1, species2, species3, k -> 4 indices
        // quaternary (4): species1, species2, species3, species4, k -> 5 indices
        int nIndices = paramType + 1;

        // Read indices line
        auto indices = readTokens(file);
        if (static_cast<int>(indices.size()) < nIndices) break;

        // Parse species indices (1-based in file)
        std::vector<int> speciesIdx(paramType);
        for (int s = 0; s < paramType; ++s) {
            speciesIdx[s] = std::stoi(indices[s]);
        }

        // Last index is k (number of coefficient sets)
        int k = std::stoi(indices[paramType]);
        k = std::max(1, k);

        // Read k sets of 2 coefficients (magnetic params use 2, not 6)
        for (int v = 0; v < k; ++v) {
            auto values = readDoubles(file, 2);
            if (values.size() < 2) break;

            int iParam = p.nMagParamCS;
            if (iParam < p.iMagneticParamCS.rows()) {
                // Column 0: number of species in parameter (2, 3, or 4)
                p.iMagneticParamCS(iParam, 0) = paramType;

                // Store species indices (keep 1-based as Fortran does)
                for (int s = 0; s < paramType; ++s) {
                    p.iMagneticParamCS(iParam, 1 + s) = speciesIdx[s];
                }

                // Set the exponent (0, 1, 2, ... for each coefficient set)
                // For ternary/quaternary, Fortran sets this to speciesIdx[v] but simplified here
                p.iMagneticParamCS(iParam, 1 + paramType) = v;

                // Store magnetic coefficients (only 2 values)
                p.dMagneticParamCS(iParam, 0) = values[0];
                p.dMagneticParamCS(iParam, 1) = values[1];
            }

            ++p.nMagParamCS;
        }
    }

    p.nMagParamPhaseCS(phaseIndex + 1) = p.nMagParamCS;

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
        // Strip Windows line endings
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        std::istringstream iss(line);
        std::string token;
        while (iss >> token) {
            // Strip any remaining \r from token
            if (!token.empty() && token.back() == '\r') {
                token.pop_back();
            }
            if (!token.empty()) {
                tokens.push_back(token);
            }
        }
    }

    return tokens;
}

std::vector<std::string> ChemSageParser::readTokens(std::ifstream& file, int count) {
    std::vector<std::string> tokens;
    std::string line;

    while (static_cast<int>(tokens.size()) < count && std::getline(file, line)) {
        // Strip Windows line endings
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        std::istringstream iss(line);
        std::string token;
        while (iss >> token && static_cast<int>(tokens.size()) < count) {
            // Strip any remaining \r from token
            if (!token.empty() && token.back() == '\r') {
                token.pop_back();
            }
            if (!token.empty()) {
                tokens.push_back(token);
            }
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
