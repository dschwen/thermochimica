/// @file QKTO.cpp
/// @brief Quasi-chemical Kohler-Toop (QKTO) model implementation
/// @author Converted from Fortran CompExcessGibbsEnergyQKTO.f90

#include "thermochimica/ThermoContext.hpp"
#include <cmath>
#include <vector>

namespace Thermochimica {

/// @brief Compute partial molar excess Gibbs energy for a QKTO sub-system
/// @details Implements the polynomial regular solution with Kohler interpolation.
///
/// For binary sub-system:
/// g^ex_lambda,z = L_z * x1^a * x2^b
///
/// With Kohler interpolation:
/// g^ex_lambda = (x1 + x2)^2 * g^ex_lambda,z
///
/// @param ctx The thermochimica context
/// @param iSolnIndex Solution phase index
/// @param iParam Parameter index
static void polyRegularQKTO(ThermoContext& ctx, int iSolnIndex, int iParam) {
    auto& thermo = *ctx.thermo;
    auto& gem = *ctx.gem;

    int iFirst = thermo.nSpeciesPhase[iSolnIndex - 1] + 1;
    int iLast = thermo.nSpeciesPhase[iSolnIndex];
    int nSpecies = iLast - iFirst + 1;

    // Initialize local variables
    std::vector<double> y(nSpecies + 1, 0.0);  // 1-based indexing
    std::vector<double> dPartialGParam(nSpecies + 1, 0.0);
    std::vector<bool> lAsymmetric1(nSpecies + 1, false);
    std::vector<bool> lAsymmetric2(nSpecies + 1, false);

    double xT = 0.0;
    double dXi1 = 0.0;
    double dXi2 = 0.0;

    // Get number of species in parameter and constituent indices
    int nSpeciesParam = thermo.iRegularParam(iParam, 1);
    int a = thermo.iRegularParam(iParam, 2);
    int b = thermo.iRegularParam(iParam, 3);
    int p = thermo.iRegularParam(iParam, nSpeciesParam + 2);
    int q = thermo.iRegularParam(iParam, nSpeciesParam + 3);

    // Set ternary parameters initially to 0
    int c = 0;
    int r = 0;
    if (nSpeciesParam == 3) {
        c = thermo.iRegularParam(iParam, 4);
        r = thermo.iRegularParam(iParam, nSpeciesParam + 4);
    }

    // Get chemical groups of first two species
    // a, b are now 0-based within phase; nSpeciesPhase[iSolnIndex - 1] is cumulative count
    int globalA = thermo.nSpeciesPhase[iSolnIndex - 1] + a;
    int globalB = thermo.nSpeciesPhase[iSolnIndex - 1] + b;

    int iGroup1 = static_cast<int>(thermo.dQKTOParams(globalA, 2));
    int iGroup2 = static_cast<int>(thermo.dQKTOParams(globalB, 2));

    // Calculate symmetry of all species with these first two
    // Note: a, b, c are 0-based, but local arrays y, lAsymmetric use 1-based indexing
    int a1 = a + 1;  // Convert to 1-based for local array access
    int b1 = b + 1;
    int c1 = c + 1;  // c is 0 if binary, so c1 would be 1 (unused for binary)

    lAsymmetric1[a1] = true;
    lAsymmetric2[b1] = true;

    // Check symmetry and handle interpolation overrides
    for (int j = thermo.nSpeciesPhase[iSolnIndex - 1] + 1;
         j <= thermo.nSpeciesPhase[iSolnIndex]; ++j) {
        int i = j - thermo.nSpeciesPhase[iSolnIndex - 1];
        int j0 = j - 1;  // 0-based index for array access

        // Check if this ternary is an exception
        bool lIsException = false;
        for (int k = 1; k <= thermo.nInterpolationOverride[iSolnIndex]; ++k) {
            bool allMatch = true;
            for (int l = 1; l <= 3; ++l) {
                int overrideIdx = thermo.iInterpolationOverride[iSolnIndex](k, l);
                if (!(overrideIdx == a || overrideIdx == b || overrideIdx == i)) {
                    allMatch = false;
                    break;
                }
            }
            if (allMatch) {
                lIsException = true;
                if (thermo.iInterpolationOverride[iSolnIndex](k, 5) == b) {
                    lAsymmetric1[i] = true;
                }
                if (thermo.iInterpolationOverride[iSolnIndex](k, 5) == a) {
                    lAsymmetric2[i] = true;
                }
                break;
            }
        }

        if (lIsException) continue;

        // If groups are unequal, assign asymmetry based on chemical group
        if (iGroup1 != iGroup2) {
            int iGroupTemp = static_cast<int>(thermo.dQKTOParams(j0, 2));
            if (iGroupTemp == iGroup1) {
                lAsymmetric1[i] = true;
            } else if (iGroupTemp == iGroup2) {
                lAsymmetric2[i] = true;
            }
        }
    }

    // Compute sum of equivalent fractions
    // Note: j is 1-based loop counter (Fortran style), use j-1 for 0-based array access
    for (int j = thermo.nSpeciesPhase[iSolnIndex - 1] + 1;
         j <= thermo.nSpeciesPhase[iSolnIndex]; ++j) {
        int j0 = j - 1;  // 0-based index for array access
        xT += thermo.dMolFraction[j0] * thermo.dQKTOParams(j0, 1);
    }

    if (xT <= 0.0) return;

    // Compute equivalent mole fractions
    for (int j = thermo.nSpeciesPhase[iSolnIndex - 1] + 1;
         j <= thermo.nSpeciesPhase[iSolnIndex]; ++j) {
        int i = j - thermo.nSpeciesPhase[iSolnIndex - 1];
        int j0 = j - 1;  // 0-based index for array access
        y[i] = thermo.dMolFraction[j0] * thermo.dQKTOParams(j0, 1) / xT;
    }

    // Calculate xis
    for (int j = thermo.nSpeciesPhase[iSolnIndex - 1] + 1;
         j <= thermo.nSpeciesPhase[iSolnIndex]; ++j) {
        int i = j - thermo.nSpeciesPhase[iSolnIndex - 1];
        if (lAsymmetric1[i]) {
            dXi1 += y[i];
        } else if (lAsymmetric2[i]) {
            dXi2 += y[i];
        }
    }

    double dXiDen = dXi1 + dXi2;

    // For ternary parameters, enforce symmetric treatment
    if (nSpeciesParam == 3) {
        c1 = c + 1;  // Update c1 for ternary case
        dXi1 = y[a1];
        dXi2 = y[b1];
        dXiDen = y[a1] + y[b1] + y[c1];
        // Force reset of all symmetry arrays
        std::fill(lAsymmetric1.begin(), lAsymmetric1.end(), false);
        std::fill(lAsymmetric2.begin(), lAsymmetric2.end(), false);
        lAsymmetric1[a1] = true;
        lAsymmetric2[b1] = true;
    }

    // Calculate g^excess for binary part
    double dGex = thermo.dExcessGibbsParam[iParam] *
                  std::pow(dXi1, p - 1) * std::pow(dXi2, q - 1) /
                  std::pow(dXiDen, p + q + r - nSpeciesParam);
    dGex *= y[a1] * y[b1] * xT;

    // Include ternary factor if present
    if (nSpeciesParam == 3) {
        dGex *= std::pow(y[c1], r - 1);
        dGex *= y[c1];
    }

    // Compute partial derivatives
    for (int j = thermo.nSpeciesPhase[iSolnIndex - 1] + 1;
         j <= thermo.nSpeciesPhase[iSolnIndex]; ++j) {
        int i = j - thermo.nSpeciesPhase[iSolnIndex - 1];
        int j0 = j - 1;  // 0-based index for array access
        double dDgex = -1.0 / xT;

        if (lAsymmetric1[i]) {
            dDgex += (p - 1) / dXi1 - (p + q + r - nSpeciesParam) / dXiDen;
        } else if (lAsymmetric2[i]) {
            dDgex += (q - 1) / dXi2 - (p + q + r - nSpeciesParam) / dXiDen;
        }

        if (i == a1) {
            dDgex += 1.0 / y[a1] / xT;
        } else if (i == b1) {
            dDgex += 1.0 / y[b1] / xT;
        } else if (i == c1 && nSpeciesParam == 3) {
            dDgex += 1.0 / y[c1] / xT;
        }

        // Ternary part of derivative
        if (nSpeciesParam == 3) {
            dDgex -= 1.0 / xT;
            if (i == c1) {
                dDgex += (r - 1) / y[c1] - (p + q + r - nSpeciesParam) / dXiDen;
            }
        }

        gem.dPartialExcessGibbs[j0] += dDgex * dGex * thermo.dQKTOParams(j0, 1);
    }
}

/// @brief Compute partial molar excess Gibbs energy for QKTO phases
/// @details Implements the Quasi-chemical Kohler-Toop model for non-ideal
/// solution phases. Uses asymmetric Kohler interpolation based on chemical groups.
///
/// @param ctx The thermochimica context
/// @param iSolnIndex Absolute index of the solution phase
void compExcessGibbsEnergyQKTO(ThermoContext& ctx, int iSolnIndex) {
    auto& thermo = *ctx.thermo;

    // Return if no interaction parameters for this phase or wrong type
    // Note: iSolnIndex is 1-based, but cSolnPhaseType is 0-based
    if ((thermo.nParamPhase[iSolnIndex] - thermo.nParamPhase[iSolnIndex - 1] == 0) ||
        (thermo.cSolnPhaseType[iSolnIndex - 1] != "QKTO")) {
        return;
    }

    // Loop through all interaction parameters in this phase
    // Note: nParamPhase is cumulative, use 0-based indexing for iParam
    int iParamStart = (iSolnIndex > 1) ? thermo.nParamPhase[iSolnIndex - 1] : 0;
    int iParamEnd = thermo.nParamPhase[iSolnIndex];

    for (int iParam = iParamStart; iParam < iParamEnd; ++iParam) {
        // Compute partial molar excess Gibbs energy for each sub-system
        // Pass iParam as 0-based index
        polyRegularQKTO(ctx, iSolnIndex, iParam);
    }
}

/// @brief Wrapper for backwards compatibility
void computeExcessGibbsQKTO(ThermoContext& ctx, int phaseIndex) {
    compExcessGibbsEnergyQKTO(ctx, phaseIndex);
}

} // namespace Thermochimica
