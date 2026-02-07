/// @file IParser.hpp
/// @brief Interface for thermodynamic database file parsers
/// @details Defines the contract for parsing different database formats

#pragma once

#include <string>
#include <vector>

namespace Thermochimica {

// Forward declarations
class ThermoState;
class ThermoIO;
class ParserState;

/// @brief Abstract interface for database file parsers
/// @details Implementations parse different database formats:
/// - ChemSageParser: ChemSage .dat format (current implementation)
/// - TDBParser: TDB format (future)
/// - JSONParser: JSON format (future)
class IParser {
public:
    virtual ~IParser() = default;

    /// @brief Parse database file into thermodynamic state
    /// @param filename Path to database file
    /// @param state Thermodynamic state to populate
    /// @param io Input/output parameters
    /// @param parserState Temporary parser workspace
    /// @return Error code (0 = success)
    virtual int parse(const std::string& filename,
                     ThermoState& state,
                     ThermoIO& io,
                     ParserState& parserState) = 0;

    /// @brief Get supported file extensions
    /// @return Vector of file extensions (e.g., {".dat", ".cs"})
    virtual std::vector<std::string> getSupportedExtensions() const = 0;

    /// @brief Check if parser can handle this file
    /// @param filename File path to check
    /// @return true if parser supports this file format
    virtual bool canParse(const std::string& filename) const = 0;

    /// @brief Get parser name for logging
    /// @return Parser name string
    virtual const char* getParserName() const = 0;
};

} // namespace Thermochimica
