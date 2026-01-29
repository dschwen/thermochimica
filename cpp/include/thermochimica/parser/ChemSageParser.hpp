#pragma once

#include "../ThermoContext.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

namespace Thermochimica {

/// ChemSage data file parser
/// Converts from ModuleParseCS.f90 pattern to C++
class ChemSageParser {
public:
    /// Parse a ChemSage data file
    /// @param ctx The thermochimica context to populate
    /// @param filename Path to the .dat file
    /// @return Error code (0 = success)
    static int parse(ThermoContext& ctx, const std::string& filename);

    /// Parse using filename already set in context
    static int parse(ThermoContext& ctx);

private:
    /// Parse the header section of the data file
    static int parseHeader(ThermoContext& ctx, std::ifstream& file);

    /// Parse the data block section
    static int parseDataBlock(ThermoContext& ctx, std::ifstream& file);

    /// Parse Gibbs energy coefficients for a species
    static int parseGibbsEquations(ThermoContext& ctx, std::ifstream& file,
                                   int speciesIndex, int nGibbsEqs, int iGibbsEqType);

    /// Parse solution phase data
    /// @param phaseType Phase type (QKTO phases need K/p parameters)
    /// @param phaseIndex Phase index (0-based) for SUBQ/SUBG coefficient storage
    /// @param localSpeciesIndex Local species index within phase (0-based)
    static int parseSolutionPhase(ThermoContext& ctx, std::ifstream& file,
                                  int speciesIndex, Constants::PhaseType phaseType,
                                  int phaseIndex = -1, int localSpeciesIndex = 0);

    /// Parse SUBL/CEF phase specific data
    /// @param hasPreHeader Output: true if stoichiometry header was found before species
    static int parseSUBLPhase(ThermoContext& ctx, std::ifstream& file,
                              int phaseIndex, bool& hasPreHeader);

    /// Parse SUBG/MQM phase specific data (before species)
    static int parseSUBGPhase(ThermoContext& ctx, std::ifstream& file,
                              int phaseIndex);

    /// Parse SUBG/SUBQ sublattice structure and mixing parameters (after species)
    static int parseSUBGExcessData(ThermoContext& ctx, std::ifstream& file,
                                   int phaseIndex, Constants::PhaseType phaseType);

    /// Parse SUBI (ionic sublattice) excess data (comes after species)
    static int parseSUBIExcessData(ThermoContext& ctx, std::ifstream& file,
                                   int phaseIndex, Constants::PhaseType phaseType);

    /// Parse mixing parameters (loop until terminator)
    static int parseMixingParametersLoop(ThermoContext& ctx, std::ifstream& file,
                                         int phaseIndex, Constants::PhaseType phaseType);

    /// Parse mixing parameters (fixed count - legacy)
    static int parseMixingParameters(ThermoContext& ctx, std::ifstream& file,
                                     int phaseIndex, int nParams);

    /// Parse magnetic parameters
    static int parseMagneticParameters(ThermoContext& ctx, std::ifstream& file,
                                       int phaseIndex);

    /// Read a line and parse tokens
    static std::vector<std::string> readTokens(std::ifstream& file);

    /// Read tokens across multiple lines until count is reached
    static std::vector<std::string> readTokens(std::ifstream& file, int count);

    /// Read a line and parse as doubles
    static std::vector<double> readDoubles(std::ifstream& file, int count);

    /// Read a line and parse as integers
    static std::vector<int> readInts(std::ifstream& file, int count);

    /// Trim whitespace from string
    static std::string trim(const std::string& str);

    /// Transfer parsed data from parser state to main state
    static void transferToThermoState(ThermoContext& ctx);
};

} // namespace Thermochimica
