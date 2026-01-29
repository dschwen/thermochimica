/// @file HeatCapacity.cpp
/// @brief Heat capacity, entropy, and enthalpy calculations

#include "thermochimica/Thermochimica.hpp"
#include "thermochimica/context/ParserState.hpp"
#include <cmath>

namespace Thermochimica {

/// @brief Compute standard heat capacity for a single species from Gibbs coefficients
/// @details For G = A1 + A2*T + A3*T*ln(T) + A4*T² + A5*T³ + A6/T + additional terms
/// Heat capacity: Cp = -T * d²G/dT² = -A3 - 2*A4*T - 6*A5*T² - 2*A6/T²
static double computeSpeciesCp(const ThermoContext& ctx, int speciesIdx) {
    auto& parser = *ctx.parser;
    auto& io = *ctx.io;

    double T = io.dTemperature;

    // Get number of Gibbs equations for this species
    int nEqs = (speciesIdx < parser.nGibbsEqSpecies.size()) ?
               parser.nGibbsEqSpecies(speciesIdx) : 1;
    if (nEqs <= 0) nEqs = 1;

    // Find the correct temperature range equation
    int iEq = 0;
    for (int eq = 0; eq < nEqs - 1; ++eq) {
        int baseIdx = eq * ParserState::kNumGibbsCoeff;
        double Tmax = parser.dGibbsCoeffSpeciesTemp(speciesIdx, baseIdx);
        if (T <= Tmax) {
            iEq = eq;
            break;
        }
        iEq = eq + 1;
    }

    // Get base index for this equation's coefficients
    int baseIdx = iEq * ParserState::kNumGibbsCoeff;

    // Get coefficients A1-A6 (stored at indices 1-6, index 0 is Tmax)
    double A3 = parser.dGibbsCoeffSpeciesTemp(speciesIdx, baseIdx + 3);  // T*ln(T) term
    double A4 = parser.dGibbsCoeffSpeciesTemp(speciesIdx, baseIdx + 4);  // T² term
    double A5 = parser.dGibbsCoeffSpeciesTemp(speciesIdx, baseIdx + 5);  // T³ term
    double A6 = parser.dGibbsCoeffSpeciesTemp(speciesIdx, baseIdx + 6);  // 1/T term

    // Compute Cp from main polynomial terms
    // Cp = -A3 - 2*A4*T - 6*A5*T² - 2*A6/T²
    double Cp = -A3 - 2.0 * A4 * T - 6.0 * A5 * T * T - 2.0 * A6 / (T * T);

    // Add contributions from additional power terms (coeff * T^exp)
    // For d²G/dT² of coeff * T^exp: coeff * exp * (exp-1) * T^(exp-2)
    // Cp contribution: -T * coeff * exp * (exp-1) * T^(exp-2) = -coeff * exp * (exp-1) * T^(exp-1)
    for (int addTerm = 0; addTerm < 3; ++addTerm) {
        int coeffIdx = baseIdx + 7 + addTerm * 2;
        int expIdx = coeffIdx + 1;

        if (expIdx < parser.dGibbsCoeffSpeciesTemp.cols()) {
            double coeff = parser.dGibbsCoeffSpeciesTemp(speciesIdx, coeffIdx);
            double exp = parser.dGibbsCoeffSpeciesTemp(speciesIdx, expIdx);

            if (std::abs(coeff) > 1e-30) {
                if (std::abs(exp - 99.0) < 0.5) {
                    // Exponent 99 means ln(T) term: coeff * ln(T)
                    // d²G/dT² = -coeff/T²
                    // Cp contribution = -T * (-coeff/T²) = coeff/T
                    Cp += coeff / T;
                } else if (std::abs(exp) > 1e-10) {
                    // General power term
                    // Cp contribution = -coeff * exp * (exp-1) * T^(exp-1)
                    Cp += -coeff * exp * (exp - 1.0) * std::pow(T, exp - 1.0);
                }
            }
        }
    }

    return Cp;
}

/// @brief Compute standard entropy for a single species from Gibbs coefficients
/// @details S = -dG/dT = -A2 - A3*(1 + ln(T)) - 2*A4*T - 3*A5*T² + A6/T²
static double computeSpeciesEntropy(const ThermoContext& ctx, int speciesIdx) {
    auto& parser = *ctx.parser;
    auto& io = *ctx.io;

    double T = io.dTemperature;
    double logT = std::log(T);

    // Get number of Gibbs equations for this species
    int nEqs = (speciesIdx < parser.nGibbsEqSpecies.size()) ?
               parser.nGibbsEqSpecies(speciesIdx) : 1;
    if (nEqs <= 0) nEqs = 1;

    // Find the correct temperature range equation
    int iEq = 0;
    for (int eq = 0; eq < nEqs - 1; ++eq) {
        int baseIdx = eq * ParserState::kNumGibbsCoeff;
        double Tmax = parser.dGibbsCoeffSpeciesTemp(speciesIdx, baseIdx);
        if (T <= Tmax) {
            iEq = eq;
            break;
        }
        iEq = eq + 1;
    }

    // Get base index for this equation's coefficients
    int baseIdx = iEq * ParserState::kNumGibbsCoeff;

    // Get coefficients A1-A6 (stored at indices 1-6, index 0 is Tmax)
    double A2 = parser.dGibbsCoeffSpeciesTemp(speciesIdx, baseIdx + 2);  // T term
    double A3 = parser.dGibbsCoeffSpeciesTemp(speciesIdx, baseIdx + 3);  // T*ln(T) term
    double A4 = parser.dGibbsCoeffSpeciesTemp(speciesIdx, baseIdx + 4);  // T² term
    double A5 = parser.dGibbsCoeffSpeciesTemp(speciesIdx, baseIdx + 5);  // T³ term
    double A6 = parser.dGibbsCoeffSpeciesTemp(speciesIdx, baseIdx + 6);  // 1/T term

    // Compute S from main polynomial terms
    // S = -dG/dT = -A2 - A3*(1 + ln(T)) - 2*A4*T - 3*A5*T² + A6/T²
    double S = -A2 - A3 * (1.0 + logT) - 2.0 * A4 * T - 3.0 * A5 * T * T + A6 / (T * T);

    // Add contributions from additional power terms
    for (int addTerm = 0; addTerm < 3; ++addTerm) {
        int coeffIdx = baseIdx + 7 + addTerm * 2;
        int expIdx = coeffIdx + 1;

        if (expIdx < parser.dGibbsCoeffSpeciesTemp.cols()) {
            double coeff = parser.dGibbsCoeffSpeciesTemp(speciesIdx, coeffIdx);
            double exp = parser.dGibbsCoeffSpeciesTemp(speciesIdx, expIdx);

            if (std::abs(coeff) > 1e-30) {
                if (std::abs(exp - 99.0) < 0.5) {
                    // ln(T) term: dG/dT = coeff/T, so S contribution = -coeff/T
                    S += -coeff / T;
                } else if (std::abs(exp) > 1e-10) {
                    // General power term: dG/dT = coeff * exp * T^(exp-1)
                    // S contribution = -coeff * exp * T^(exp-1)
                    S += -coeff * exp * std::pow(T, exp - 1.0);
                }
            }
        }
    }

    return S;
}

void computeHeatCapacity(ThermoContext& ctx) {
    auto& io = *ctx.io;
    auto& thermo = *ctx.thermo;

    if (!io.lHeatCapacityEntropyEnthalpy) {
        return;
    }

    double T = io.dTemperature;
    double R = Constants::kIdealGasConstant;

    // Compute system heat capacity as sum over all species
    // Cp_sys = Σ (n_i * Cp_i)
    double totalCp = 0.0;
    double totalS = 0.0;
    double totalH = 0.0;

    for (int i = 0; i < thermo.nSpecies; ++i) {
        double moles = thermo.dMolesSpecies(i);
        if (moles > 1e-100) {
            // Standard heat capacity for this species
            double Cp_i = computeSpeciesCp(ctx, i);
            totalCp += moles * Cp_i;

            // Standard entropy for this species
            double S_i = computeSpeciesEntropy(ctx, i);

            // For solution phase species, add mixing entropy contribution
            // S_mix = -R * Σ x_i * ln(x_i)
            if (thermo.dMolFraction(i) > 1e-100 && thermo.dMolFraction(i) < 1.0 - 1e-10) {
                S_i -= R * std::log(thermo.dMolFraction(i));
            }
            totalS += moles * S_i;

            // Enthalpy: H = G + T*S
            double G_i = thermo.dStdGibbsEnergy(i);
            double H_i = G_i + T * S_i;
            totalH += moles * H_i;
        }
    }

    io.dHeatCapacity = totalCp;
    io.dEntropy = totalS;
    io.dEnthalpy = totalH;
}

} // namespace Thermochimica
