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
    static int parseSolutionPhase(ThermoContext& ctx, std::ifstream& file,
                                  int phaseIndex);

    /// Parse SUBL/CEF phase specific data
    static int parseSUBLPhase(ThermoContext& ctx, std::ifstream& file,
                              int phaseIndex);

    /// Parse SUBG/MQM phase specific data
    static int parseSUBGPhase(ThermoContext& ctx, std::ifstream& file,
                              int phaseIndex);

    /// Parse mixing parameters
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
