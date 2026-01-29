/// @file RKMP.cpp
/// @brief Redlich-Kister-Muggiano-Polynomial (RKMP/RKMPM) model implementation
/// @author Converted from Fortran CompExcessGibbsEnergyRKMP.f90

#include "thermochimica/ThermoContext.hpp"
#include <cmath>
#include <algorithm>

namespace Thermochimica {

/// @brief Compute partial molar excess Gibbs energy for RKMP/RKMPM phases
/// @details Implements the Redlich-Kister-Muggianu polynomial model for
/// non-ideal solution phases with binary, ternary, and quaternary interactions.
///
/// For binary interactions (index z):
/// g^ex_lambda,z = x1 * x2 * sum_{v=0}^{Nz} (^v)L_{1,2} * (x1 - x2)^v
///
/// For ternary interactions:
/// g^ex_lambda,z = x1 * x2 * x3 * ((1-x1-x2-x3)/3 + xj) * L_{1,2,3}
///
/// For quaternary interactions:
/// g^ex_lambda,z = x1 * x2 * x3 * x4 * L_{1,2,3,4}
///
/// @param ctx The thermochimica context
/// @param iSolnIndex Absolute index of the solution phase
void compExcessGibbsEnergyRKMP(ThermoContext& ctx, int iSolnIndex) {
    auto& thermo = *ctx.thermo;
    auto& gem = *ctx.gem;
    auto& io = *ctx.io;

    // Return if no interaction parameters for this phase
    if (thermo.nParamPhase[iSolnIndex] - thermo.nParamPhase[iSolnIndex - 1] == 0) {
        return;
    }

    // Only proceed if the phase type is correct
    if (iSolnIndex < 0 || iSolnIndex >= static_cast<int>(thermo.iSolnPhaseType.size())) {
        return;
    }
    Constants::PhaseType phaseType = thermo.iSolnPhaseType[iSolnIndex];
    if (phaseType != Constants::PhaseType::RKMP && phaseType != Constants::PhaseType::RKMPM) {
        return;
    }

    // Store indices of first and last species in this phase (1-based like Fortran)
    int iFirstSpecies = thermo.nSpeciesPhase[iSolnIndex - 1] + 1;
    int iLastSpecies = thermo.nSpeciesPhase[iSolnIndex];

    // Loop through all interaction parameters in this phase
    for (int iParam = thermo.nParamPhase[iSolnIndex - 1] + 1;
         iParam <= thermo.nParamPhase[iSolnIndex]; ++iParam) {

        // Compute temporary variables for convenience
        int idx1 = iFirstSpecies + thermo.iRegularParam(iParam, 2) - 1;
        int idx2 = iFirstSpecies + thermo.iRegularParam(iParam, 3) - 1;

        double x1 = thermo.dMolFraction[idx1];
        double x2 = thermo.dMolFraction[idx2];
        double xprod = x1 * x2;
        double dx = x1 - x2;

        int paramType = thermo.iRegularParam(iParam, 1);

        if (paramType == 2) {
            // ======================
            // Binary parameter
            // ======================

            // Skip if dx = 0 to prevent INF or NAN
            if (dx == 0.0) continue;

            int iExponent = thermo.iRegularParam(iParam, 4);
            double dxvmo = std::pow(dx, iExponent - 1);
            dxvmo = std::min(dxvmo, 1.0e30);  // Prevent overflow

            // Loop through species in this phase
            for (int i = iFirstSpecies; i <= iLastSpecies; ++i) {
                int j = i - iFirstSpecies + 1;

                if (j == thermo.iRegularParam(iParam, 2)) {
                    // First species of parameter
                    gem.dPartialExcessGibbs[i] += thermo.dExcessGibbsParam[iParam] *
                        dxvmo * ((x2 - xprod) * dx +
                                 static_cast<double>(iExponent) * xprod * (1.0 - dx));
                }
                else if (j == thermo.iRegularParam(iParam, 3)) {
                    // Second species of parameter
                    gem.dPartialExcessGibbs[i] += thermo.dExcessGibbsParam[iParam] *
                        dxvmo * ((x1 - xprod) * dx -
                                 static_cast<double>(iExponent) * xprod * (1.0 + dx));
                }
                else {
                    // This species does not belong to the parameter
                    gem.dPartialExcessGibbs[i] -= thermo.dExcessGibbsParam[iParam] *
                        std::pow(dx, iExponent) * xprod *
                        (1.0 + static_cast<double>(iExponent));
                }
            }
        }
        else if (paramType == 3) {
            // ======================
            // Ternary parameter
            // ======================
            int idx3 = iFirstSpecies + thermo.iRegularParam(iParam, 4) - 1;
            double x3 = thermo.dMolFraction[idx3];
            int jIdx = thermo.iRegularParam(iParam, 5) + iFirstSpecies - 1;
            double xj = thermo.dMolFraction[jIdx];
            xprod = xprod * x3;

            // Loop through species in this phase
            for (int i = iFirstSpecies; i <= iLastSpecies; ++i) {
                int j = i - iFirstSpecies + 1;

                // Compute Kronecker-Delta term for ternary parameter
                double KD = 0.0;
                if (thermo.iRegularParam(iParam, 5) == j) {
                    KD = 1.0;
                }

                if (j == thermo.iRegularParam(iParam, 2) ||
                    j == thermo.iRegularParam(iParam, 3) ||
                    j == thermo.iRegularParam(iParam, 4)) {
                    // This species contributes to the parameter
                    double xi = thermo.dMolFraction[i];
                    gem.dPartialExcessGibbs[i] += thermo.dExcessGibbsParam[iParam] *
                        xprod * (((1.0 / xi) - 3.0) *
                                 ((1.0 - x1 - x2 - x3) / 3.0 + xj) + KD);
                }
                else {
                    // This species does not contribute to the parameter
                    gem.dPartialExcessGibbs[i] += thermo.dExcessGibbsParam[iParam] *
                        xprod * (x1 + x2 + x3 - 3.0 * xj - 2.0 / 3.0);
                }
            }
        }
        else if (paramType == 4) {
            // ======================
            // Quaternary parameter
            // ======================
            int idx3 = iFirstSpecies + thermo.iRegularParam(iParam, 4) - 1;
            int idx4 = iFirstSpecies + thermo.iRegularParam(iParam, 5) - 1;
            double x3 = thermo.dMolFraction[idx3];
            double x4 = thermo.dMolFraction[idx4];
            xprod = xprod * x3 * x4;

            // Loop through species in this phase
            for (int i = iFirstSpecies; i <= iLastSpecies; ++i) {
                int j = i - iFirstSpecies + 1;

                if (j == thermo.iRegularParam(iParam, 2) ||
                    j == thermo.iRegularParam(iParam, 3) ||
                    j == thermo.iRegularParam(iParam, 4) ||
                    j == thermo.iRegularParam(iParam, 5)) {
                    // This species contributes to the parameter
                    double xi = thermo.dMolFraction[i];
                    gem.dPartialExcessGibbs[i] += thermo.dExcessGibbsParam[iParam] *
                        xprod * ((1.0 / xi) - 3.0);
                }
                else {
                    // This species does not contribute to the parameter
                    gem.dPartialExcessGibbs[i] -= 3.0 * xprod *
                        thermo.dExcessGibbsParam[iParam];
                }
            }
        }
        else {
            // Parameter index not supported/recognized
            io.INFOThermo = 32;
            return;
        }
    }
}

/// @brief Wrapper for backwards compatibility
void computeExcessGibbsRKMP(ThermoContext& ctx, int phaseIndex) {
    compExcessGibbsEnergyRKMP(ctx, phaseIndex);
}

} // namespace Thermochimica
