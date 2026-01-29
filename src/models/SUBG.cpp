/// @file SUBG.cpp
/// @brief Modified Quasichemical Model (SUBG/SUBQ) implementation
/// @author Converted from Fortran CompExcessGibbsEnergySUBG.f90

#include "thermochimica/ThermoContext.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

namespace Thermochimica {

/// @brief Compute chemical potentials for SUBG/SUBQ phases (Modified Quasichemical Model)
/// @details Implements the Modified Quasichemical Model for phases with short-range order.
/// The model focuses on pairs of nearest neighbors (quadruplets) rather than species.
///
/// References:
/// - A.D. Pelton et al., Met. Trans. B, 31B (2000) 651-659
/// - A.D. Pelton, P. Chartrand, Met. Trans. B, 32A (2001) 1355-1360
///
/// @param ctx The thermochimica context
/// @param iSolnIndex Absolute index of the solution phase
void compExcessGibbsEnergySUBG(ThermoContext& ctx, int iSolnIndex) {
    auto& thermo = *ctx.thermo;
    auto& gem = *ctx.gem;
    auto& io = *ctx.io;

    // Only proceed for SUBG or SUBQ phases
    if (iSolnIndex < 0 || iSolnIndex >= static_cast<int>(thermo.iSolnPhaseType.size())) {
        return;
    }
    Constants::PhaseType phaseType = thermo.iSolnPhaseType[iSolnIndex];
    if (phaseType != Constants::PhaseType::SUBG && phaseType != Constants::PhaseType::SUBQ) {
        return;
    }

    // Define temporary variables
    int iFirst = thermo.nSpeciesPhase[iSolnIndex - 1] + 1;
    int iLast = thermo.nSpeciesPhase[iSolnIndex];
    int iSPI = thermo.iPhaseSublattice[iSolnIndex];
    int nSub1 = thermo.nConstituentSublattice(iSPI, 1);
    int nSub2 = thermo.nConstituentSublattice(iSPI, 2);
    int nA2X2 = nSub1 * nSub2;
    int nPhaseElements = nSub1 + nSub2;

    // Allocate arrays
    std::vector<double> dXi(nPhaseElements + 1, 0.0);
    std::vector<double> dYi(nPhaseElements + 1, 0.0);
    std::vector<double> dFi(nPhaseElements + 1, 0.0);
    std::vector<double> dNi(nPhaseElements + 1, 0.0);
    std::vector<std::vector<double>> dXij(nSub1 + 1, std::vector<double>(nSub2 + 1, 0.0));
    std::vector<std::vector<double>> dNij(nSub1 + 1, std::vector<double>(nSub2 + 1, 0.0));
    std::vector<std::vector<double>> dXsij(nSub1 + 1, std::vector<double>(nSub2 + 1, 0.0));
    std::vector<std::vector<double>> dNsij(nSub1 + 1, std::vector<double>(nSub2 + 1, 0.0));
    std::vector<bool> lAsymmetric1(std::max(nSub1, nSub2) + 1, false);
    std::vector<bool> lAsymmetric2(std::max(nSub1, nSub2) + 1, false);

    // Initialize chemical potentials
    for (int i = iFirst; i <= iLast; ++i) {
        thermo.dChemicalPotential[i] = 0.0;
        gem.dPartialExcessGibbs[i] = 0.0;
    }

    // ==================================
    // Compute X_i and Y_i
    // ==================================

    // Do cations first
    double dSum = 0.0;
    for (int i = 1; i <= nSub1; ++i) {
        for (int k = 1; k <= thermo.nPairsSRO(iSPI, 2); ++k) {
            int l = iFirst + k - 1;
            double dZa = thermo.dCoordinationNumber[iSPI]( k, 1);
            double dZb = thermo.dCoordinationNumber[iSPI]( k, 2);
            if (i == thermo.iPairID[iSPI]( k, 1)) {
                dNi[i] += thermo.dMolFraction[l] / dZa;
                dYi[i] += thermo.dMolFraction[l] / 2.0;
            }
            if (i == thermo.iPairID[iSPI]( k, 2)) {
                dNi[i] += thermo.dMolFraction[l] / dZb;
                dYi[i] += thermo.dMolFraction[l] / 2.0;
            }
        }
        dSum += dNi[i];
    }
    for (int i = 1; i <= nSub1; ++i) {
        dXi[i] = dNi[i] / dSum;
    }

    // Do anions
    dSum = 0.0;
    for (int i = 1; i <= nSub2; ++i) {
        int j = i + nSub1;
        for (int k = 1; k <= thermo.nPairsSRO(iSPI, 2); ++k) {
            int l = iFirst + k - 1;
            double dZx = thermo.dCoordinationNumber[iSPI]( k, 3);
            double dZy = thermo.dCoordinationNumber[iSPI]( k, 4);
            if (j == thermo.iPairID[iSPI]( k, 3)) {
                dNi[j] += thermo.dMolFraction[l] / dZx;
                dYi[j] += thermo.dMolFraction[l] / 2.0;
            }
            if (j == thermo.iPairID[iSPI]( k, 4)) {
                dNi[j] += thermo.dMolFraction[l] / dZy;
                dYi[j] += thermo.dMolFraction[l] / 2.0;
            }
        }
        dSum += dNi[j];
    }
    for (int i = 1; i <= nSub2; ++i) {
        int j = i + nSub1;
        dXi[j] = dNi[j] / dSum;
    }

    // ==================================
    // Compute X_i/j
    // ==================================
    double dSumNij = 0.0;
    double dSumNsij = 0.0;
    for (int i = 1; i <= nSub1; ++i) {
        for (int j = 1; j <= nSub2; ++j) {
            for (int k = 1; k <= thermo.nPairsSRO(iSPI, 2); ++k) {
                int l = iFirst + k - 1;
                int nA = 0;
                if (i == thermo.iPairID[iSPI]( k, 1)) nA++;
                if (i == thermo.iPairID[iSPI]( k, 2)) nA++;
                int nX = 0;
                if ((j + nSub1) == thermo.iPairID[iSPI]( k, 3)) nX++;
                if ((j + nSub1) == thermo.iPairID[iSPI]( k, 4)) nX++;

                int m = 1;
                for (int kk = 1; kk <= nA2X2; ++kk) {
                    if (thermo.iConstituentSublattice[iSPI]( 1, kk) == i &&
                        thermo.iConstituentSublattice[iSPI]( 2, kk) == j) {
                        m = kk;
                    }
                }
                dNij[i][j] += thermo.dMolFraction[l] * nA * nX;
                dNsij[i][j] += thermo.dMolFraction[l] * nA * nX / thermo.dZetaSpecies(iSPI, m);
            }
            dSumNij += dNij[i][j];
            dSumNsij += dNsij[i][j];
        }
    }

    for (int i = 1; i <= nSub1; ++i) {
        for (int j = 1; j <= nSub2; ++j) {
            dXij[i][j] = dNij[i][j] / dSumNij;
            dXsij[i][j] = dNsij[i][j] / dSumNsij;
        }
    }

    // Calculate F_i
    for (int i = 1; i <= nSub1; ++i) {
        for (int j = 1; j <= nSub2; ++j) {
            dFi[i] += dXsij[i][j];
            int k = j + nSub1;
            dFi[k] += dXsij[i][j];
        }
    }

    // ==================================
    // COMPUTE REFERENCE GIBBS ENERGY AND IDEAL MIXING TERMS
    // ==================================
    for (int k = 1; k <= thermo.nPairsSRO(iSPI, 2); ++k) {
        double dConfEntropy = 0.0;
        int l = iFirst + k - 1;

        // Coordination numbers for this quadruplet
        double dZa = thermo.dCoordinationNumber[iSPI]( k, 1);
        double dZb = thermo.dCoordinationNumber[iSPI]( k, 2);
        double dZx = thermo.dCoordinationNumber[iSPI]( k, 3);
        double dZy = thermo.dCoordinationNumber[iSPI]( k, 4);

        // Loop over n_i contributions to entropy - Cations
        for (int i = 1; i <= nSub1; ++i) {
            if (i == thermo.iPairID[iSPI]( k, 1) && dXi[i] > 0.0) {
                dConfEntropy += std::log(dXi[i]) / dZa;
            }
            if (i == thermo.iPairID[iSPI]( k, 2) && dXi[i] > 0.0) {
                dConfEntropy += std::log(dXi[i]) / dZb;
            }
        }
        // Anions
        for (int i = 1; i <= nSub2; ++i) {
            int j = i + nSub1;
            if (j == thermo.iPairID[iSPI]( k, 3) && dXi[j] > 0.0) {
                dConfEntropy += std::log(dXi[j]) / dZx;
            }
            if (j == thermo.iPairID[iSPI]( k, 4) && dXi[j] > 0.0) {
                dConfEntropy += std::log(dXi[j]) / dZy;
            }
        }

        // Loop over n_i/j contributions to entropy
        for (int i = 1; i <= nSub1; ++i) {
            for (int j = 1; j <= nSub2; ++j) {
                int nA = 0;
                if (i == thermo.iPairID[iSPI]( k, 1)) nA++;
                if (i == thermo.iPairID[iSPI]( k, 2)) nA++;
                int nX = 0;
                if ((j + nSub1) == thermo.iPairID[iSPI]( k, 3)) nX++;
                if ((j + nSub1) == thermo.iPairID[iSPI]( k, 4)) nX++;

                int m = 1;
                for (int kk = 1; kk <= nA2X2; ++kk) {
                    if (thermo.iConstituentSublattice[iSPI]( 1, kk) == i &&
                        thermo.iConstituentSublattice[iSPI]( 2, kk) == j) {
                        m = kk;
                    }
                }

                double denom = dFi[i] * dFi[j + nSub1];
                if (dXsij[i][j] > 0.0 && denom > 0.0) {
                    dConfEntropy += std::log(dXsij[i][j] / denom) *
                                   (nA * nX / thermo.dZetaSpecies(iSPI, m));
                }
            }
        }

        // Pair indices
        int ii = thermo.iPairID[iSPI]( k, 1);
        int jj = thermo.iPairID[iSPI]( k, 2);
        int kk = thermo.iPairID[iSPI]( k, 3);
        int ll = thermo.iPairID[iSPI]( k, 4);
        int ka = kk - nSub1;
        int la = ll - nSub1;

        // Add n_ij/kl contribution
        int iWeight = 1;
        if (ii != jj) iWeight *= 2;
        if (kk != ll) iWeight *= 2;

        // SUBG and SUBQ differ in entropy calculation
        double dPowXij = 0.0, dPowYi = 0.0;
        if (phaseType == Constants::PhaseType::SUBG) {
            dPowXij = 1.0;
            dPowYi = 1.0;
        } else if (phaseType == Constants::PhaseType::SUBQ) {
            dPowXij = 0.75;
            dPowYi = 0.5;
        }

        double yProd = dYi[ii] * dYi[jj] * dYi[kk] * dYi[ll];
        if (yProd != 0.0) {
            dSum = iWeight *
                   std::pow(dXij[ii][ka] * dXij[ii][la] * dXij[jj][ka] * dXij[jj][la], dPowXij) /
                   std::pow(yProd, dPowYi);
            if (dSum == 0.0) {
                dConfEntropy = 100.0;
            } else if (thermo.dMolFraction[l] > 0.0) {
                dConfEntropy += std::log(thermo.dMolFraction[l] / dSum);
            }
        }

        double dRef = thermo.dStdGibbsEnergy[l];
        thermo.dChemicalPotential[l] = dRef + dConfEntropy;
    }

    // ==================================
    // EXCESS MIXING PARAMETERS
    // ==================================
    for (int abxy = thermo.nParamPhase[iSolnIndex - 1] + 1;
         abxy <= thermo.nParamPhase[iSolnIndex]; ++abxy) {

        if (thermo.dExcessGibbsParam[abxy] == 0.0) continue;

        // AB/XY parametrization
        int a = thermo.iRegularParam(abxy, 2);
        int b = thermo.iRegularParam(abxy, 3);
        int xx = thermo.iRegularParam(abxy, 4);
        int yy = thermo.iRegularParam(abxy, 5);
        int x = xx - nSub1;
        int y = yy - nSub1;
        double p = thermo.iRegularParam(abxy, 6);
        double q = thermo.iRegularParam(abxy, 7);
        // Note: r, s, d, w from iRegularParam indices 8-11 are reserved for future use

        // Calculate block index
        int iBlock;
        if (x == y) {
            iBlock = (x - 1) * (nSub1 * (nSub1 + 1) / 2);
        } else if (x > y) {
            continue;  // Skip
        } else {
            iBlock = (nSub2 + (x - 1) + ((y - 2) * (y - 1) / 2)) *
                     (nSub1 * (nSub1 + 1) / 2);
        }
        if (a == b) {
            iBlock = iBlock + a;
        } else if (a > b) {
            continue;  // Skip
        } else {
            iBlock = iBlock + nSub1 + a + ((b - 2) * (b - 1) / 2);
        }
        iBlock = iBlock + iFirst - 1;

        // Calculate xi and chi for symmetry
        double dXi1 = 0.0, dXi2 = 0.0, dXiDen = 0.0;
        double dChi1 = 0.0, dChi2 = 0.0, dChiDen = 0.0;
        std::fill(lAsymmetric1.begin(), lAsymmetric1.end(), false);
        std::fill(lAsymmetric2.begin(), lAsymmetric2.end(), false);

        // Setup asymmetry based on binary type
        if (x == y) {
            // Cation-cation mixing with common anion
            lAsymmetric1[a] = true;
            lAsymmetric2[b] = true;

            // Check chemical groups for asymmetry
            for (int i = 1; i <= nSub1; ++i) {
                if (thermo.iChemicalGroup[iSPI]( 1, a) != thermo.iChemicalGroup[iSPI]( 1, b)) {
                    if (thermo.iChemicalGroup[iSPI]( 1, i) == thermo.iChemicalGroup[iSPI]( 1, a)) {
                        lAsymmetric1[i] = true;
                    } else if (thermo.iChemicalGroup[iSPI]( 1, i) == thermo.iChemicalGroup[iSPI]( 1, b)) {
                        lAsymmetric2[i] = true;
                    }
                }
            }

            // Calculate chi terms
            for (int ijkl = 1; ijkl <= thermo.nPairsSRO(iSPI, 2); ++ijkl) {
                int i = thermo.iPairID[iSPI]( ijkl, 1);
                int j = thermo.iPairID[iSPI]( ijkl, 2);
                int k_idx = thermo.iPairID[iSPI]( ijkl, 3) - nSub1;
                int l_idx = thermo.iPairID[iSPI]( ijkl, 4) - nSub1;
                int iQuad = ijkl + iFirst - 1;

                if (lAsymmetric1[i] && lAsymmetric1[j]) {
                    if (x == k_idx && x == l_idx) {
                        dChi1 += thermo.dMolFraction[iQuad];
                    } else if (phaseType == Constants::PhaseType::SUBQ && (x == k_idx || x == l_idx)) {
                        dChi1 += 0.5 * thermo.dMolFraction[iQuad];
                    }
                }
                if (lAsymmetric2[i] && lAsymmetric2[j]) {
                    if (x == k_idx && x == l_idx) {
                        dChi2 += thermo.dMolFraction[iQuad];
                    } else if (phaseType == Constants::PhaseType::SUBQ && (x == k_idx || x == l_idx)) {
                        dChi2 += 0.5 * thermo.dMolFraction[iQuad];
                    }
                }
                if ((lAsymmetric1[i] || lAsymmetric2[i]) &&
                    (lAsymmetric1[j] || lAsymmetric2[j])) {
                    if (x == k_idx && x == l_idx) {
                        dChiDen += thermo.dMolFraction[iQuad];
                    } else if (phaseType == Constants::PhaseType::SUBQ && (x == k_idx || x == l_idx)) {
                        dChiDen += 0.5 * thermo.dMolFraction[iQuad];
                    }
                }
            }

            // Calculate xi terms
            for (int i = 1; i <= nSub1; ++i) {
                for (int k = 1; k <= thermo.nPairsSRO(iSPI, 2); ++k) {
                    int l = k + iFirst - 1;
                    double contrib = thermo.dMolFraction[l] / 4.0;
                    if (lAsymmetric1[i]) {
                        if (i == thermo.iPairID[iSPI]( k, 1) && xx == thermo.iPairID[iSPI]( k, 3))
                            dXi1 += contrib;
                        if (i == thermo.iPairID[iSPI]( k, 1) && xx == thermo.iPairID[iSPI]( k, 4))
                            dXi1 += contrib;
                        if (i == thermo.iPairID[iSPI]( k, 2) && xx == thermo.iPairID[iSPI]( k, 3))
                            dXi1 += contrib;
                        if (i == thermo.iPairID[iSPI]( k, 2) && xx == thermo.iPairID[iSPI]( k, 4))
                            dXi1 += contrib;
                    }
                    if (lAsymmetric2[i]) {
                        if (i == thermo.iPairID[iSPI]( k, 1) && xx == thermo.iPairID[iSPI]( k, 3))
                            dXi2 += contrib;
                        if (i == thermo.iPairID[iSPI]( k, 1) && xx == thermo.iPairID[iSPI]( k, 4))
                            dXi2 += contrib;
                        if (i == thermo.iPairID[iSPI]( k, 2) && xx == thermo.iPairID[iSPI]( k, 3))
                            dXi2 += contrib;
                        if (i == thermo.iPairID[iSPI]( k, 2) && xx == thermo.iPairID[iSPI]( k, 4))
                            dXi2 += contrib;
                    }
                }
            }
        } else if (a == b) {
            // Anion-anion mixing with common cation
            lAsymmetric1[x] = true;
            lAsymmetric2[y] = true;

            if (thermo.iChemicalGroup[iSPI]( 2, x) != thermo.iChemicalGroup[iSPI]( 2, y)) {
                for (int i = 1; i <= nSub2; ++i) {
                    if (thermo.iChemicalGroup[iSPI]( 2, i) == thermo.iChemicalGroup[iSPI]( 2, x)) {
                        lAsymmetric1[i] = true;
                    } else if (thermo.iChemicalGroup[iSPI]( 2, i) == thermo.iChemicalGroup[iSPI]( 2, y)) {
                        lAsymmetric2[i] = true;
                    }
                }
            }

            // Calculate chi and xi similar to above for anion mixing
            for (int ijkl = 1; ijkl <= thermo.nPairsSRO(iSPI, 2); ++ijkl) {
                int i = thermo.iPairID[iSPI]( ijkl, 1);
                int j = thermo.iPairID[iSPI]( ijkl, 2);
                int k_idx = thermo.iPairID[iSPI]( ijkl, 3) - nSub1;
                int l_idx = thermo.iPairID[iSPI]( ijkl, 4) - nSub1;
                int iQuad = ijkl + iFirst - 1;

                if (lAsymmetric1[k_idx] && lAsymmetric1[l_idx]) {
                    if (a == i && a == j) {
                        dChi1 += thermo.dMolFraction[iQuad];
                    } else if (phaseType == Constants::PhaseType::SUBQ && (a == i || a == j)) {
                        dChi1 += 0.5 * thermo.dMolFraction[iQuad];
                    }
                }
                if (lAsymmetric2[k_idx] && lAsymmetric2[l_idx]) {
                    if (a == i && a == j) {
                        dChi2 += thermo.dMolFraction[iQuad];
                    } else if (phaseType == Constants::PhaseType::SUBQ && (a == i || a == j)) {
                        dChi2 += 0.5 * thermo.dMolFraction[iQuad];
                    }
                }
                if ((lAsymmetric1[k_idx] || lAsymmetric2[k_idx]) &&
                    (lAsymmetric1[l_idx] || lAsymmetric2[l_idx])) {
                    if (a == i && a == j) {
                        dChiDen += thermo.dMolFraction[iQuad];
                    } else if (phaseType == Constants::PhaseType::SUBQ && (a == i || a == j)) {
                        dChiDen += 0.5 * thermo.dMolFraction[iQuad];
                    }
                }
            }

            for (int i = 1; i <= nSub2; ++i) {
                int ii = i + nSub1;
                for (int k = 1; k <= thermo.nPairsSRO(iSPI, 2); ++k) {
                    int l = k + iFirst - 1;
                    double contrib = thermo.dMolFraction[l] / 4.0;
                    if (lAsymmetric1[i]) {
                        if (a == thermo.iPairID[iSPI]( k, 1) && ii == thermo.iPairID[iSPI]( k, 3))
                            dXi1 += contrib;
                        if (a == thermo.iPairID[iSPI]( k, 1) && ii == thermo.iPairID[iSPI]( k, 4))
                            dXi1 += contrib;
                        if (a == thermo.iPairID[iSPI]( k, 2) && ii == thermo.iPairID[iSPI]( k, 3))
                            dXi1 += contrib;
                        if (a == thermo.iPairID[iSPI]( k, 2) && ii == thermo.iPairID[iSPI]( k, 4))
                            dXi1 += contrib;
                    }
                    if (lAsymmetric2[i]) {
                        if (a == thermo.iPairID[iSPI]( k, 1) && ii == thermo.iPairID[iSPI]( k, 3))
                            dXi2 += contrib;
                        if (a == thermo.iPairID[iSPI]( k, 1) && ii == thermo.iPairID[iSPI]( k, 4))
                            dXi2 += contrib;
                        if (a == thermo.iPairID[iSPI]( k, 2) && ii == thermo.iPairID[iSPI]( k, 3))
                            dXi2 += contrib;
                        if (a == thermo.iPairID[iSPI]( k, 2) && ii == thermo.iPairID[iSPI]( k, 4))
                            dXi2 += contrib;
                    }
                }
            }
        }

        dXiDen = dXi1 + dXi2;
        if (dChiDen > 0.0) {
            dChi1 /= dChiDen;
            dChi2 /= dChiDen;
        }

        // Calculate excess energy contribution
        double dGex = 0.0;
        double dDgexBase = 0.0;
        char paramType = thermo.cRegularParam[abxy];

        // G-type terms
        if (paramType == 'G') {
            dGex = thermo.dExcessGibbsParam[abxy] * std::pow(dChi1, p) * std::pow(dChi2, q);
            dDgexBase = -dGex * (p + q) / dChiDen;
        }
        // Q-type terms
        else if (paramType == 'Q') {
            if (dXiDen > 0.0) {
                dGex = thermo.dExcessGibbsParam[abxy] *
                       std::pow(dXi1, p) * std::pow(dXi2, q) /
                       std::pow(dXiDen, p + q);
                dDgexBase = -dGex * (p + q) / dXiDen;
            }
        }
        // B-type terms
        else if (paramType == 'B' || paramType == 'H') {
            double dXtot = dXsij[a][x] + dXsij[b][y];
            if (dXtot > 0.0) {
                dGex = thermo.dExcessGibbsParam[abxy] *
                       std::pow(dXsij[a][x], 1.0 + p) *
                       std::pow(dXsij[b][y], 1.0 + q) /
                       std::pow(dXtot, 1.0 + p + q);
            }
            // B-type has different derivative handling
            // Simplified for now
        }
        // R-type (reciprocal) terms
        else if (paramType == 'R') {
            gem.dPartialExcessGibbs[iBlock] += dGex;
            continue;
        }
        else {
            io.INFOThermo = 42;
            continue;
        }

        // Add contribution to block quadruplet
        gem.dPartialExcessGibbs[iBlock] += dGex / 2.0;

        // Add contributions from related quadruplets
        if (a == b && x != y) {
            // Add g^ex to quads AC/XY
            for (int c = 1; c <= nSub1; ++c) {
                if (c == a) continue;
                int e = std::min(a, c);
                int f = std::max(a, c);
                int ia = (a < c) ? 1 : 2;
                int iQuad = (nSub2 + (x - 1) + ((y - 2) * (y - 1) / 2)) *
                           (nSub1 * (nSub1 + 1) / 2) +
                           nSub1 + e + ((f - 2) * (f - 1) / 2);
                iQuad += iFirst - 1;
                double ratio = thermo.dCoordinationNumber[iSPI]( iBlock - iFirst + 1, 1) /
                              thermo.dCoordinationNumber[iSPI]( iQuad - iFirst + 1, ia);
                gem.dPartialExcessGibbs[iQuad] += (dGex / 4.0) * ratio;
            }
        }

        if (a != b && x == y) {
            // Add g^ex to quads AB/XZ
            for (int z = 1; z <= nSub2; ++z) {
                if (z == x) continue;
                int e = std::min(x, z);
                int f = std::max(x, z);
                int ix = (x < z) ? 3 : 4;
                int iQuad = (nSub2 + (e - 1) + ((f - 2) * (f - 1) / 2)) *
                           (nSub1 * (nSub1 + 1) / 2) +
                           nSub1 + a + ((b - 2) * (b - 1) / 2);
                iQuad += iFirst - 1;
                double ratio = thermo.dCoordinationNumber[iSPI]( iBlock - iFirst + 1, 3) /
                              thermo.dCoordinationNumber[iSPI]( iQuad - iFirst + 1, ix);
                gem.dPartialExcessGibbs[iQuad] += (dGex / 4.0) * ratio;
            }
        }

        // Add dg^ex contributions to all quadruplets
        for (int ijkl = 1; ijkl <= thermo.nPairsSRO(iSPI, 2); ++ijkl) {
            int iQuad2 = ijkl + iFirst - 1;
            int i = thermo.iPairID[iSPI]( ijkl, 1);
            int j = thermo.iPairID[iSPI]( ijkl, 2);
            int k_idx = thermo.iPairID[iSPI]( ijkl, 3) - nSub1;
            int l_idx = thermo.iPairID[iSPI]( ijkl, 4) - nSub1;

            double dChiFactor = 0.0;
            double dDgex = 0.0;

            if (a != b && x == y) {
                if (x == k_idx && x == l_idx) {
                    dChiFactor = 1.0;
                } else if (phaseType == Constants::PhaseType::SUBQ && (x == k_idx || x == l_idx)) {
                    dChiFactor = 0.5;
                }
            } else if (a == b && x != y) {
                if (a == i && a == j) {
                    dChiFactor = 1.0;
                } else if (phaseType == Constants::PhaseType::SUBQ && (a == i || a == j)) {
                    dChiFactor = 0.5;
                }
            }

            // Calculate derivatives based on parameter type
            if (paramType == 'G') {
                if (lAsymmetric1[i] && lAsymmetric1[j] && dChi1 > 0.0) {
                    dDgex += dChiFactor * dGex * p / dChi1 / dChiDen;
                }
                if (lAsymmetric2[i] && lAsymmetric2[j] && dChi2 > 0.0) {
                    dDgex += dChiFactor * dGex * q / dChi2 / dChiDen;
                }
                if ((lAsymmetric1[i] || lAsymmetric2[i]) &&
                    (lAsymmetric1[j] || lAsymmetric2[j])) {
                    dDgex += dChiFactor * dDgexBase;
                }
            } else if (paramType == 'Q') {
                for (int ii = 1; ii <= nSub1; ++ii) {
                    if (lAsymmetric1[ii]) {
                        if (ii == i && x == k_idx)
                            dDgex += dDgexBase / 4 + dGex * p / (4 * dXi1);
                        if (ii == i && x == l_idx)
                            dDgex += dDgexBase / 4 + dGex * p / (4 * dXi1);
                        if (ii == j && x == k_idx)
                            dDgex += dDgexBase / 4 + dGex * p / (4 * dXi1);
                        if (ii == j && x == l_idx)
                            dDgex += dDgexBase / 4 + dGex * p / (4 * dXi1);
                    }
                    if (lAsymmetric2[ii]) {
                        if (ii == i && x == k_idx)
                            dDgex += dDgexBase / 4 + dGex * q / (4 * dXi2);
                        if (ii == i && x == l_idx)
                            dDgex += dDgexBase / 4 + dGex * q / (4 * dXi2);
                        if (ii == j && x == k_idx)
                            dDgex += dDgexBase / 4 + dGex * q / (4 * dXi2);
                        if (ii == j && x == l_idx)
                            dDgex += dDgexBase / 4 + dGex * q / (4 * dXi2);
                    }
                }
            }

            gem.dPartialExcessGibbs[iQuad2] += thermo.dMolFraction[iBlock] * dDgex / 2.0;
        }
    }
}

/// @brief Wrapper for backwards compatibility
void computeExcessGibbsSUBG(ThermoContext& ctx, int phaseIndex) {
    compExcessGibbsEnergySUBG(ctx, phaseIndex);
}

} // namespace Thermochimica
