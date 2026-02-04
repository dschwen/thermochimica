/// @file batch_validate.cpp
/// @brief Batch validation tool for Thermochimica
/// @details Reads JSON test cases, runs calculations, outputs results for validation
///
/// Usage: batch_validate <input.json> [output.json] [--compare reference.json]

#include "thermochimica/ThermoClass.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <map>

using namespace Thermochimica;

/// @brief Simple JSON parser (minimal implementation)
/// @details This is a minimal JSON implementation to avoid external dependencies.
/// For production use, consider nlohmann/json or similar library.
struct JsonValue {
    enum Type { OBJECT, ARRAY, STRING, NUMBER, BOOLEAN, NULL_TYPE };
    Type type;
    std::string stringValue;
    double numberValue;
    bool boolValue;
    std::map<std::string, JsonValue> objectValue;
    std::vector<JsonValue> arrayValue;

    JsonValue() : type(NULL_TYPE), numberValue(0.0), boolValue(false) {}

    static JsonValue parseNumber(const std::string& str) {
        JsonValue v;
        v.type = NUMBER;
        v.numberValue = std::stod(str);
        return v;
    }

    static JsonValue parseString(const std::string& str) {
        JsonValue v;
        v.type = STRING;
        v.stringValue = str;
        return v;
    }

    static JsonValue parseBoolean(bool b) {
        JsonValue v;
        v.type = BOOLEAN;
        v.boolValue = b;
        return v;
    }

    static JsonValue parseNull() {
        JsonValue v;
        v.type = NULL_TYPE;
        return v;
    }
};

/// @brief Test case structure
struct TestCase {
    std::string name;
    std::string database;
    double temperature;
    double pressure;
    std::map<std::string, double> composition;

    // Optional reference values
    bool hasReference = false;
    double refGibbs = 0.0;
    std::map<std::string, double> refPhases;
    std::map<std::string, double> refChemPot;
};

/// @brief Test result structure
struct TestResult {
    std::string name;
    bool success;
    int errorCode;
    std::string errorMessage;
    double gibbs;
    std::map<std::string, double> phases;
    std::map<std::string, double> chemPot;

    // Comparison metrics (if reference provided)
    bool compared = false;
    double gibbsError = 0.0;
    std::map<std::string, double> phaseErrors;
    std::map<std::string, double> chemPotErrors;
};

/// @brief Parse element composition from JSON-like string
/// @details Expected format: {"Element1": value1, "Element2": value2, ...}
std::map<std::string, double> parseComposition(const std::string& str) {
    std::map<std::string, double> comp;
    size_t pos = str.find('{');
    if (pos == std::string::npos) return comp;

    std::string content = str.substr(pos + 1);
    pos = content.find('}');
    if (pos != std::string::npos) {
        content = content.substr(0, pos);
    }

    // Parse "key": value pairs
    size_t start = 0;
    while (start < content.length()) {
        // Find key
        size_t keyStart = content.find('"', start);
        if (keyStart == std::string::npos) break;
        size_t keyEnd = content.find('"', keyStart + 1);
        if (keyEnd == std::string::npos) break;

        std::string key = content.substr(keyStart + 1, keyEnd - keyStart - 1);

        // Find colon
        size_t colon = content.find(':', keyEnd);
        if (colon == std::string::npos) break;

        // Find value
        size_t valStart = colon + 1;
        while (valStart < content.length() && std::isspace(content[valStart])) valStart++;

        size_t valEnd = content.find_first_of(",}", valStart);
        if (valEnd == std::string::npos) valEnd = content.length();

        std::string valStr = content.substr(valStart, valEnd - valStart);
        double value = std::stod(valStr);

        comp[key] = value;
        start = valEnd + 1;
    }

    return comp;
}

/// @brief Read test cases from JSON file
std::vector<TestCase> readTestCases(const std::string& filename) {
    std::vector<TestCase> cases;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open input file: " << filename << std::endl;
        return cases;
    }

    std::string line;
    TestCase currentCase;
    bool inTestCase = false;

    while (std::getline(file, line)) {
        // Trim whitespace
        size_t start = line.find_first_not_of(" \t\r\n");
        size_t end = line.find_last_not_of(" \t\r\n");
        if (start == std::string::npos) continue;
        line = line.substr(start, end - start + 1);

        if (line.find("\"name\":") != std::string::npos) {
            if (inTestCase && !currentCase.name.empty()) {
                cases.push_back(currentCase);
            }
            currentCase = TestCase();
            inTestCase = true;

            size_t pos = line.find(':');
            if (pos != std::string::npos) {
                std::string val = line.substr(pos + 1);
                size_t q1 = val.find('"');
                size_t q2 = val.find('"', q1 + 1);
                if (q1 != std::string::npos && q2 != std::string::npos) {
                    currentCase.name = val.substr(q1 + 1, q2 - q1 - 1);
                }
            }
        }
        else if (line.find("\"database\":") != std::string::npos) {
            size_t pos = line.find(':');
            if (pos != std::string::npos) {
                std::string val = line.substr(pos + 1);
                size_t q1 = val.find('"');
                size_t q2 = val.find('"', q1 + 1);
                if (q1 != std::string::npos && q2 != std::string::npos) {
                    currentCase.database = val.substr(q1 + 1, q2 - q1 - 1);
                }
            }
        }
        else if (line.find("\"temperature\":") != std::string::npos) {
            size_t pos = line.find(':');
            if (pos != std::string::npos) {
                std::string val = line.substr(pos + 1);
                size_t comma = val.find(',');
                if (comma != std::string::npos) val = val.substr(0, comma);
                currentCase.temperature = std::stod(val);
            }
        }
        else if (line.find("\"pressure\":") != std::string::npos) {
            size_t pos = line.find(':');
            if (pos != std::string::npos) {
                std::string val = line.substr(pos + 1);
                size_t comma = val.find(',');
                if (comma != std::string::npos) val = val.substr(0, comma);
                currentCase.pressure = std::stod(val);
            }
        }
        else if (line.find("\"composition\":") != std::string::npos) {
            // Read composition object (may span multiple lines)
            std::string compStr = line;
            while (compStr.find('}') == std::string::npos) {
                if (!std::getline(file, line)) break;
                compStr += line;
            }
            currentCase.composition = parseComposition(compStr);
        }
        else if (line.find("\"reference\":") != std::string::npos) {
            currentCase.hasReference = true;
            // Read reference object
            while (std::getline(file, line)) {
                if (line.find("\"gibbs\":") != std::string::npos) {
                    size_t pos = line.find(':');
                    if (pos != std::string::npos) {
                        std::string val = line.substr(pos + 1);
                        size_t comma = val.find(',');
                        if (comma != std::string::npos) val = val.substr(0, comma);
                        currentCase.refGibbs = std::stod(val);
                    }
                }
                if (line.find('}') != std::string::npos) break;
            }
        }
    }

    // Add last case
    if (inTestCase && !currentCase.name.empty()) {
        cases.push_back(currentCase);
    }

    file.close();
    return cases;
}

/// @brief Run a single test case
TestResult runTestCase(const TestCase& testCase) {
    TestResult result;
    result.name = testCase.name;

    try {
        ThermoClass thermo;

        // Load database
        int loadResult = thermo.loadDatabase(testCase.database);
        if (loadResult != 0) {
            result.success = false;
            result.errorCode = loadResult;
            result.errorMessage = thermo.getErrorMessage();
            return result;
        }

        // Set conditions
        thermo.setStandardUnits();
        thermo.setTemperaturePressure(testCase.temperature, testCase.pressure);

        for (const auto& [element, mass] : testCase.composition) {
            thermo.setElementMass(element, mass);
        }

        // Run calculation
        int calcResult = thermo.calculate();
        result.errorCode = calcResult;

        if (calcResult != 0) {
            result.success = false;
            result.errorMessage = thermo.getErrorMessage();
            return result;
        }

        result.success = true;
        result.gibbs = thermo.getGibbsEnergy();

        // Get phase amounts (we don't know phase names ahead of time, so skip for now)
        // In a real implementation, you'd query the database for phase names

        // Get chemical potentials for elements in composition
        for (const auto& [element, mass] : testCase.composition) {
            auto [chemPot, info] = thermo.getOutputChemPot(element);
            if (info == 0) {
                result.chemPot[element] = chemPot;
            }
        }

        // Compare with reference if available
        if (testCase.hasReference) {
            result.compared = true;
            result.gibbsError = std::abs(result.gibbs - testCase.refGibbs);
        }

    } catch (const std::exception& e) {
        result.success = false;
        result.errorCode = -999;
        result.errorMessage = std::string("Exception: ") + e.what();
    }

    return result;
}

/// @brief Write results to JSON file
void writeResults(const std::string& filename, const std::vector<TestResult>& results) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open output file: " << filename << std::endl;
        return;
    }

    file << std::setprecision(16);
    file << "{\n";
    file << "  \"results\": [\n";

    for (size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];

        file << "    {\n";
        file << "      \"name\": \"" << r.name << "\",\n";
        file << "      \"success\": " << (r.success ? "true" : "false") << ",\n";
        file << "      \"error_code\": " << r.errorCode << ",\n";

        if (!r.success) {
            file << "      \"error_message\": \"" << r.errorMessage << "\",\n";
        }

        file << "      \"gibbs_energy\": " << r.gibbs << ",\n";

        // Chemical potentials
        file << "      \"chemical_potentials\": {\n";
        size_t cpCount = 0;
        for (const auto& [elem, val] : r.chemPot) {
            file << "        \"" << elem << "\": " << val;
            if (++cpCount < r.chemPot.size()) file << ",";
            file << "\n";
        }
        file << "      }";

        // Comparison metrics
        if (r.compared) {
            file << ",\n";
            file << "      \"comparison\": {\n";
            file << "        \"gibbs_error\": " << r.gibbsError << "\n";
            file << "      }";
        }

        file << "\n    }";
        if (i < results.size() - 1) file << ",";
        file << "\n";
    }

    file << "  ]\n";
    file << "}\n";

    file.close();
}

/// @brief Print summary to console
void printSummary(const std::vector<TestResult>& results) {
    int passed = 0;
    int failed = 0;
    double maxGibbsError = 0.0;
    std::string maxErrorCase;

    for (const auto& r : results) {
        if (r.success) {
            passed++;
            if (r.compared && r.gibbsError > maxGibbsError) {
                maxGibbsError = r.gibbsError;
                maxErrorCase = r.name;
            }
        } else {
            failed++;
        }
    }

    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Total tests: " << results.size() << std::endl;
    std::cout << "Passed: " << passed << std::endl;
    std::cout << "Failed: " << failed << std::endl;

    if (maxGibbsError > 0.0) {
        std::cout << "\nComparison results:" << std::endl;
        std::cout << "Max Gibbs error: " << maxGibbsError << " J" << std::endl;
        std::cout << "  in test case: " << maxErrorCase << std::endl;
    }

    if (failed > 0) {
        std::cout << "\nFailed tests:" << std::endl;
        for (const auto& r : results) {
            if (!r.success) {
                std::cout << "  " << r.name << ": " << r.errorMessage << std::endl;
            }
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input.json> [output.json]" << std::endl;
        std::cerr << "\nReads test cases from input JSON file and runs Thermochimica calculations." << std::endl;
        std::cerr << "Results are written to output JSON file (default: results.json)" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = (argc >= 3) ? argv[2] : "results.json";

    std::cout << "Reading test cases from: " << inputFile << std::endl;

    auto testCases = readTestCases(inputFile);

    if (testCases.empty()) {
        std::cerr << "Error: No test cases found in input file" << std::endl;
        return 1;
    }

    std::cout << "Found " << testCases.size() << " test cases" << std::endl;
    std::cout << "\nRunning calculations..." << std::endl;

    std::vector<TestResult> results;

    for (size_t i = 0; i < testCases.size(); ++i) {
        const auto& tc = testCases[i];
        std::cout << "  [" << (i+1) << "/" << testCases.size() << "] "
                  << tc.name << "... " << std::flush;

        auto result = runTestCase(tc);
        results.push_back(result);

        if (result.success) {
            std::cout << "OK (G = " << std::scientific << std::setprecision(6)
                      << result.gibbs << " J)" << std::endl;
        } else {
            std::cout << "FAILED (code " << result.errorCode << ")" << std::endl;
        }
    }

    std::cout << "\nWriting results to: " << outputFile << std::endl;
    writeResults(outputFile, results);

    printSummary(results);

    return 0;
}
