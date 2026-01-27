/// @file SUBL.cpp
/// @brief Compound Energy Formalism (SUBL/SUBLM) model implementation
/// @author Converted from Fortran CompExcessGibbsEnergySUBL.f90

#include "thermochimica/ThermoContext.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

namespace Thermochimica {

// Forward declaration for magnetic contribution
void compGibbsMagneticSoln(ThermoContext& ctx, int iSolnIndex);

/// @brief Compute partial molar excess Gibbs energy for SUBL/SUBLM phases
/// @details Implements the Compound Energy Formalism (CEF) for phases with sublattices.
///
/// The chemical potential of a component in a SUBL phase is:
/// mu_i = sum_j x_j * g_j^o * (1 - Ns + sum_s delta_{i,j}/y_{j(s)})
///        + RT * sum_s a_s * ln(y_{i(s)}) + g_i^ex
///
/// The partial molar excess Gibbs energy is:
/// g_i^ex = sum_p (prod_m y_{m(s)}) * sum_z (^z)L_{j,k} * (1 - (Ns + z) + sum_s delta_{i,p}/y_{i(s)})
///
/// @param ctx The thermochimica context
/// @param iSolnIndex Absolute index of the solution phase
void compExcessGibbsEnergySUBL(ThermoContext& ctx, int iSolnIndex) {
    auto& thermo = *ctx.thermo;
    auto& gem = *ctx.gem;
    auto& io = *ctx.io;

    // Only proceed if the correct phase type is selected
    const auto& phaseType = thermo.cSolnPhaseType[iSolnIndex];
    if (phaseType != "SUBL" && phaseType != "SUBLM") {
        return;
    }

    // Define temporary variables
    int iChargedPhaseID = thermo.iPhaseSublattice[iSolnIndex];
    int nSublattice = thermo.nSublatticePhase[iChargedPhaseID];
    int iFirst = thermo.nSpeciesPhase[iSolnIndex - 1] + 1;
    int iLast = thermo.nSpeciesPhase[iSolnIndex];

    // Initialize arrays
    // Reset site fractions
    for (int s = 1; s <= nSublattice; ++s) {
        for (int c = 1; c <= thermo.nMaxConstituentSys; ++c) {
            thermo.dSiteFraction[iChargedPhaseID]( s, c) = 0.0;
        }
    }

    // Reset chemical potentials and partial excess Gibbs
    for (int i = iFirst; i <= iLast; ++i) {
        thermo.dChemicalPotential[i] = 0.0;
        gem.dPartialExcessGibbs[i] = 0.0;
    }

    // Vector to store sum of site fractions on each sublattice
    std::vector<double> dTempVec(thermo.nMaxSublatticeSys + 1, 0.0);

    // ==================================
    // Compute site fractions on each sublattice
    // ==================================
    for (int i = iFirst; i <= iLast; ++i) {
        int m = i - iFirst + 1;  // Relative component index

        for (int s = 1; s <= nSublattice; ++s) {
            int c = thermo.iConstituentSublattice[iChargedPhaseID]( s, m);
            thermo.dSiteFraction[iChargedPhaseID]( s, c) += thermo.dMolFraction[i];
        }
    }

    // Compute sum of site fractions on each sublattice and inverse
    for (int s = 1; s <= nSublattice; ++s) {
        for (int c = 1; c <= thermo.nMaxConstituentSys; ++c) {
            dTempVec[s] += thermo.dSiteFraction[iChargedPhaseID]( s, c);
        }
        dTempVec[s] = 1.0 / dTempVec[s];
    }

    // Correct mole fractions by site fractions
    for (int i = iFirst; i <= iLast; ++i) {
        double dTemp = 1.0;
        int m = i - iFirst + 1;

        for (int s = 1; s <= nSublattice; ++s) {
            int c = thermo.iConstituentSublattice[iChargedPhaseID]( s, m);
            dTemp *= thermo.dSiteFraction[iChargedPhaseID]( s, c) * dTempVec[s];
        }

        thermo.dMolFraction[i] = dTemp;
    }

    // ==================================
    // REFERENCE GIBBS ENERGY AND IDEAL MIXING
    // ==================================
    for (int i = iFirst; i <= iLast; ++i) {
        int m = i - iFirst + 1;

        // Loop through components for reference energy
        for (int j = iFirst; j <= iLast; ++j) {
            int n = j - iFirst + 1;

            // Compute pre-factor term
            double dTemp = 1.0 - static_cast<double>(nSublattice);

            // Loop through sublattices
            for (int s = 1; s <= nSublattice; ++s) {
                int k = thermo.iConstituentSublattice[iChargedPhaseID]( s, m);
                int l = thermo.iConstituentSublattice[iChargedPhaseID]( s, n);

                // Apply Kronecker-Delta term to pre-factor
                if (k == l) {
                    dTemp += 1.0 / thermo.dSiteFraction[iChargedPhaseID]( s, k);
                }
            }

            // Update reference molar Gibbs energy
            thermo.dChemicalPotential[i] += dTemp * thermo.dMolFraction[j] * thermo.dStdGibbsEnergy[j];
        }

        // Add ideal mixing contribution
        for (int s = 1; s <= nSublattice; ++s) {
            int c = thermo.iConstituentSublattice[iChargedPhaseID]( s, m);
            double siteFrac = thermo.dSiteFraction[iChargedPhaseID]( s, c);
            if (siteFrac > 0.0) {
                thermo.dChemicalPotential[i] += thermo.dStoichSublattice(iChargedPhaseID, s) *
                    std::log(siteFrac);
            }
        }

        // Sum stoichiometry and add penalty for vacancies-only species
        double dSum = 0.0;
        for (int j = 1; j <= thermo.nElements; ++j) {
            dSum += std::abs(thermo.dStoichSpecies(i, j));
        }
        if (dSum == 0.0 && thermo.dMolFraction[i] > 0.9) {
            thermo.dChemicalPotential[i] += 1000.0;
        }
    }

    // ==================================
    // MAGNETIC TERMS (for SUBLM phases)
    // ==================================
    if (phaseType == "SUBLM") {
        compGibbsMagneticSoln(ctx, iSolnIndex);
    }

    // ==================================
    // EXCESS TERMS
    // ==================================

    // Return if no interaction parameters
    if (thermo.nParamPhase[iSolnIndex] - thermo.nParamPhase[iSolnIndex - 1] == 0) {
        return;
    }

    // Loop through parameters
    for (int l = thermo.nParamPhase[iSolnIndex - 1] + 1;
         l <= thermo.nParamPhase[iSolnIndex]; ++l) {

        // Initialize temporary variables
        double dPreFactor = 1.0;
        int iFirstParam = 0, iSecondParam = 0, iThirdParam = 0;
        int iSubParam = 0, iSubParam2 = 0, iExponent = 0;
        double dFirstParam = 0.0, dSecondParam = 0.0, dThirdParam = 0.0;
        double dFirstParam2 = 0.0, dSecondParam2 = 0.0;
        int iFirstParam2 = 0, iSecondParam2 = 0;
        int iMixType = 0;

        // Store number of constituents in this parameter
        int n = thermo.iRegularParam(l, 1);

        // Determine the mixing parameter type
        if (thermo.iSUBLParamData(l, 1) == 1 && thermo.iSUBLParamData(l, 3) == 2) {
            // Binary mixing on one sublattice
            iMixType = 2;

            // Loop through constituents associated with this parameter
            for (int k = 2; k <= n + 1; ++k) {
                int c = thermo.iRegularParam(l, k) % 10000;
                int s = thermo.iRegularParam(l, k) / 10000;

                // Compute prefactor term
                dPreFactor *= thermo.dSiteFraction[iChargedPhaseID]( s, c);

                // Store first and second site fractions
                if (k == thermo.iSUBLParamData(l, 2)) {
                    dFirstParam = thermo.dSiteFraction[iChargedPhaseID]( s, c);
                    iFirstParam = c;
                    iSubParam = s;
                } else if (k == thermo.iSUBLParamData(l, 2) + 1) {
                    dSecondParam = thermo.dSiteFraction[iChargedPhaseID]( s, c);
                    iSecondParam = c;
                }
            }

            // Multiply prefactor by excess Gibbs energy parameter
            iExponent = thermo.iRegularParam(l, n + 2);
            dPreFactor *= thermo.dExcessGibbsParam[l] *
                          std::pow(dFirstParam - dSecondParam, iExponent);
        }
        else if (thermo.iSUBLParamData(l, 1) == 1 && thermo.iSUBLParamData(l, 3) == 3) {
            // Ternary mixing on one sublattice
            iMixType = 3;
            int iTernaryCon = thermo.iRegularParam(l, n + 2);

            for (int k = 2; k <= n + 1; ++k) {
                int c = thermo.iRegularParam(l, k) % 10000;
                int s = thermo.iRegularParam(l, k) / 10000;

                dPreFactor *= thermo.dSiteFraction[iChargedPhaseID]( s, c);

                if (k == thermo.iSUBLParamData(l, 2) + iTernaryCon) {
                    dFirstParam = thermo.dSiteFraction[iChargedPhaseID]( s, c);
                    iFirstParam = c;
                    iSubParam = s;
                } else if (k == thermo.iSUBLParamData(l, 2) + (iTernaryCon + 1) % 3) {
                    dSecondParam = thermo.dSiteFraction[iChargedPhaseID]( s, c);
                    iSecondParam = c;
                } else if (k == thermo.iSUBLParamData(l, 2) + (iTernaryCon + 2) % 3) {
                    dThirdParam = thermo.dSiteFraction[iChargedPhaseID]( s, c);
                    iThirdParam = c;
                }
            }

            iExponent = 0;
            dPreFactor *= thermo.dExcessGibbsParam[l];
            dPreFactor *= (dFirstParam + (1.0 - dFirstParam - dSecondParam - dThirdParam) / 3.0);
        }
        else if (thermo.iSUBLParamData(l, 1) == 2 &&
                 thermo.iSUBLParamData(l, 3) == 2 &&
                 thermo.iSUBLParamData(l, 5) == 2) {
            // Binary mixing on two sublattices (reciprocal)
            iMixType = 4;

            for (int k = 2; k <= n + 1; ++k) {
                int c = thermo.iRegularParam(l, k) % 10000;
                int s = thermo.iRegularParam(l, k) / 10000;

                dPreFactor *= thermo.dSiteFraction[iChargedPhaseID]( s, c);

                if (k == thermo.iSUBLParamData(l, 2)) {
                    dFirstParam = thermo.dSiteFraction[iChargedPhaseID]( s, c);
                    iFirstParam = c;
                    iSubParam = s;
                } else if (k == thermo.iSUBLParamData(l, 2) + 1) {
                    dSecondParam = thermo.dSiteFraction[iChargedPhaseID]( s, c);
                    iSecondParam = c;
                } else if (k == thermo.iSUBLParamData(l, 4)) {
                    dFirstParam2 = thermo.dSiteFraction[iChargedPhaseID]( s, c);
                    iFirstParam2 = c;
                    iSubParam2 = s;
                } else if (k == thermo.iSUBLParamData(l, 4) + 1) {
                    dSecondParam2 = thermo.dSiteFraction[iChargedPhaseID]( s, c);
                    iSecondParam2 = c;
                }
            }

            // Handle odd/even exponent switching
            if (thermo.iRegularParam(l, n + 2) % 2 == 0) {
                iExponent = thermo.iRegularParam(l, n + 2) / 2;
            } else {
                iExponent = (thermo.iRegularParam(l, n + 2) - 1) / 2;
                // Swap first and second parameters
                std::swap(iFirstParam, iFirstParam2);
                std::swap(iSecondParam, iSecondParam2);
                std::swap(dFirstParam, dFirstParam2);
                std::swap(dSecondParam, dSecondParam2);
                std::swap(iSubParam, iSubParam2);
            }

            dPreFactor *= thermo.dExcessGibbsParam[l] *
                          std::pow(dFirstParam - dSecondParam, iExponent);
        }
        else {
            // Unrecognized excess mixing term
            io.INFOThermo = 36;
            return;
        }

        // Loop through species in phase
        for (int i = iFirst; i <= iLast; ++i) {
            double KD = 0.0;
            int m = i - iFirst + 1;
            double dTemp;

            if (iMixType == 2 || iMixType == 4) {
                dTemp = -static_cast<double>(nSublattice + iExponent);
            } else if (iMixType == 3) {
                dTemp = -static_cast<double>(nSublattice + 1);
            }

            // Loop through sublattices
            for (int s = 1; s <= nSublattice; ++s) {
                int c = thermo.iConstituentSublattice[iChargedPhaseID]( s, m);

                // Assign Kronecker-Delta term
                if (s == iSubParam) {
                    if (c == iFirstParam) {
                        if (iMixType == 2 || iMixType == 4) {
                            KD = 1.0;
                        } else if (iMixType == 3) {
                            KD = 2.0 / 3.0;
                        }
                    } else if (c == iSecondParam) {
                        if (iMixType == 2 || iMixType == 4) {
                            KD = -1.0;
                        } else if (iMixType == 3) {
                            KD = -1.0 / 3.0;
                        }
                    } else if (c == iThirdParam) {
                        if (iMixType == 3) {
                            KD = -1.0 / 3.0;
                        }
                    } else {
                        KD = 0.0;
                    }
                    if (iMixType == 3) {
                        KD -= (dFirstParam + (-dFirstParam - dSecondParam - dThirdParam) / 3.0);
                    }
                }

                // Loop through constituents in parameter
                for (int j = 2; j <= n + 1; ++j) {
                    int k_const = thermo.iRegularParam(l, j) % 10000;
                    int d = thermo.iRegularParam(l, j) / 10000;

                    // Skip if different sublattice or constituent
                    if (d != s || k_const != c) continue;

                    // Include contribution from this site fraction
                    dTemp += 1.0 / thermo.dSiteFraction[iChargedPhaseID]( s, c);
                }
            }

            // Apply higher order terms
            if (iMixType == 2) {
                if (dFirstParam != dSecondParam) {
                    dTemp += KD * static_cast<double>(iExponent) / (dFirstParam - dSecondParam);
                }
            } else if (iMixType == 3) {
                double denom = dFirstParam + (1.0 - dFirstParam - dSecondParam - dThirdParam) / 3.0;
                if (denom != 0.0) {
                    dTemp += KD / denom;
                }
            } else if (iMixType == 4) {
                if (dFirstParam != dSecondParam) {
                    dTemp += KD * static_cast<double>(iExponent) / (dFirstParam - dSecondParam);
                }
                dTemp -= 1.0;
            }

            // Apply partial molar excess Gibbs energy
            gem.dPartialExcessGibbs[i] += dPreFactor * dTemp;
        }
    }
}

/// @brief Wrapper for backwards compatibility
void computeExcessGibbsSUBL(ThermoContext& ctx, int phaseIndex) {
    compExcessGibbsEnergySUBL(ctx, phaseIndex);
}

/// @brief Compute site fractions from species mole fractions
void computeSiteFractions(ThermoContext& ctx, int phaseIndex) {
    auto& thermo = *ctx.thermo;
    auto& gem = *ctx.gem;

    int iChargedPhaseID = thermo.iPhaseSublattice[phaseIndex];
    if (iChargedPhaseID <= 0) return;

    int nSublattice = thermo.nSublatticePhase[iChargedPhaseID];
    int iFirst = thermo.nSpeciesPhase[phaseIndex - 1] + 1;
    int iLast = thermo.nSpeciesPhase[phaseIndex];

    // Reset site fractions
    for (int s = 1; s <= nSublattice; ++s) {
        for (int c = 1; c <= thermo.nMaxConstituentSys; ++c) {
            thermo.dSiteFraction[iChargedPhaseID]( s, c) = 0.0;
        }
    }

    // Compute site fractions
    for (int i = iFirst; i <= iLast; ++i) {
        int m = i - iFirst + 1;
        for (int s = 1; s <= nSublattice; ++s) {
            int c = thermo.iConstituentSublattice[iChargedPhaseID]( s, m);
            thermo.dSiteFraction[iChargedPhaseID]( s, c) += thermo.dMolFraction[i];
        }
    }
}

/// @brief Magnetic contribution to Gibbs energy for SUBLM phases
/// @note Not yet implemented - returns zero contribution
void compGibbsMagneticSoln(ThermoContext& ctx, int iSolnIndex) {
    (void)ctx;
    (void)iSolnIndex;
    // Magnetic ordering contribution based on Curie temperature and magnetic moment
    // would be implemented here for SUBLM phases
    static bool warned = false;
    if (!warned) {
        std::cerr << "Warning: Magnetic Gibbs energy contribution (SUBLM) not yet implemented\n";
        warned = true;
    }
}

} // namespace Thermochimica
