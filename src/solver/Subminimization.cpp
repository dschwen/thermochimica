/// @file Subminimization.cpp
/// @brief Subminimization solver for non-ideal solution phase stability
/// @author Converted from Fortran Subminimization.f90

#include "thermochimica/ThermoContext.hpp"
#include "thermochimica/util/Tolerances.hpp"
#include <Eigen/LU>
#include <cmath>
#include <vector>
#include <algorithm>

namespace Thermochimica {

// Forward declarations
void compExcessGibbsEnergy(ThermoContext& ctx, int phaseIndex);

/// @brief Compute chemical potentials of solution phase constituents
/// @param iFirst 0-based first species index (inclusive)
/// @param iLast 0-based last species index (exclusive)
static void subMinChemicalPotential(ThermoContext& ctx, int iSolnPhaseIndex,
                                    int iFirst, int iLast) {
    auto& thermo = *ctx.thermo;
    auto& gem = *ctx.gem;

    // Initialize chemical potentials
    for (int i = iFirst; i < iLast; ++i) {
        thermo.dChemicalPotential[i] = 0.0;
        gem.dPartialExcessGibbs[i] = 0.0;
    }

    // Compute excess Gibbs energy contributions
    compExcessGibbsEnergy(ctx, iSolnPhaseIndex);
}

/// @brief Compute driving force for a solution phase
static double subMinDrivingForce(ThermoContext& ctx, int iFirst, int iLast,
                                 const std::vector<double>& dChemicalPotentialStar) {
    auto& thermo = *ctx.thermo;
    int nVar = iLast - iFirst + 1;
    double dDrivingForce = 0.0;

    for (int j = 0; j < nVar; ++j) {
        int i = iFirst + j;
        dDrivingForce += thermo.dMolFraction[i] *
            (thermo.dChemicalPotential[i] - dChemicalPotentialStar[j]);
    }

    return dDrivingForce / static_cast<double>(nVar);
}

/// @brief Compute functional norm for convergence testing
/// @param iFirst 0-based first species index (inclusive)
/// @param iLast 0-based last species index (exclusive)
static double subMinFunctionNorm(ThermoContext& ctx, int iSolnPhaseIndex,
                                 int iFirst, int iLast) {
    auto& thermo = *ctx.thermo;

    double dSubMinFunctionNorm = -1.0;

    // Compute residual of mole fractions
    for (int i = iFirst; i < iLast; ++i) {
        dSubMinFunctionNorm += thermo.dMolFraction[i];
    }

    // Compute residual of charge neutrality constraints if ionic
    if (thermo.iPhaseElectronID[iSolnPhaseIndex] != 0) {
        int j = thermo.iPhaseElectronID[iSolnPhaseIndex];
        for (int i = iFirst; i < iLast; ++i) {
            dSubMinFunctionNorm += thermo.dMolFraction[i] * thermo.dStoichSpecies(i, j);
        }
    }

    return std::abs(dSubMinFunctionNorm);
}

/// @brief Newton step for subminimization
static int subMinNewton(ThermoContext& ctx, int iSolnPhaseIndex,
                        int iFirst, int /*iLast*/, int nVar,
                        const std::vector<double>& dChemicalPotentialStar,
                        double dDrivingForce,
                        Eigen::VectorXd& dRHS) {
    auto& thermo = *ctx.thermo;

    // Number of equations
    int nEqn = nVar + 1;
    if (thermo.iPhaseElectronID[iSolnPhaseIndex] != 0) {
        nEqn += 1;
    }

    Eigen::MatrixXd dHessian = Eigen::MatrixXd::Zero(nEqn, nEqn);
    dRHS = Eigen::VectorXd::Constant(nEqn, -1.0);

    // Construct diagonal and part of the arrow head
    for (int j = 0; j < nVar; ++j) {
        int i = iFirst + j;
        double molFrac = thermo.dMolFraction[i];

        dHessian(j, j) = std::min(std::max(1.0 / molFrac, 1.0 + 1e-10), 1e10);
        dHessian(nVar, j) = -1.0;
        dHessian(j, nVar) = -1.0;
        dRHS(j) = dDrivingForce -
                  (thermo.dChemicalPotential[i] + 1.0 - dChemicalPotentialStar[j]);
        dRHS(nVar) += molFrac;
    }

    // Apply additional row/column if the phase is ionic
    if (thermo.iPhaseElectronID[iSolnPhaseIndex] != 0) {
        int k = thermo.iPhaseElectronID[iSolnPhaseIndex];
        int lastRow = nEqn - 1;
        dRHS(lastRow) = 0.0;

        for (int j = 0; j < nVar; ++j) {
            int i = iFirst + j;
            dHessian(lastRow, j) = -thermo.dStoichSpecies(i, k);
            dHessian(j, lastRow) = dHessian(lastRow, j);
            dRHS(lastRow) += thermo.dStoichSpecies(i, k) * thermo.dMolFraction[i];
        }
    }

    // Solve the system
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(dHessian);
    if (lu.determinant() == 0.0) {
        return -1;  // Singular matrix
    }

    dRHS = lu.solve(dRHS);
    return 0;
}

/// @brief Line search for subminimization
static void subMinLineSearch(ThermoContext& ctx, int iSolnPhaseIndex,
                             int iFirst, int iLast, int nVar,
                             const std::vector<double>& dChemicalPotentialStar,
                             Eigen::VectorXd& dRHS,
                             double& dDrivingForce,
                             double& dDrivingForceLast,
                             bool& lSubMinConverged) {
    auto& thermo = *ctx.thermo;

    const double dMaxChange = 0.2;
    const double dMinMoleFraction = 1e-100;
    const double dSubMinTolerance = 1e-10;

    double dStepLength = 1.0;

    // Initialize steplength to prevent negative mole fractions
    for (int j = 0; j < nVar; ++j) {
        int i = iFirst + j;

        if (thermo.dMolFraction[i] + dRHS(j) <= 0.0) {
            double dTemp = -0.99 * thermo.dMolFraction[i] / dRHS(j);
            dStepLength = std::min(dStepLength, dTemp);
        }

        if (dRHS(j) != 0.0) {
            double dTemp = std::abs(dMaxChange / dRHS(j));
            dStepLength = std::min(dStepLength, dTemp);
        }
    }

    // Update mole fractions
    double dMaxUpdate = 0.0;
    for (int j = 0; j < nVar; ++j) {
        int i = iFirst + j;
        double dChange = dStepLength * dRHS(j);

        if (thermo.dMolFraction[i] + dChange > 1.0) {
            dChange = 1.0 - thermo.dMolFraction[i];
        }

        thermo.dMolFraction[i] += dChange;
        dMaxUpdate = std::max(dMaxUpdate, std::abs(dChange));
    }

    // Iterate to satisfy Wolfe conditions
    for (int k = 0; k < 5; ++k) {
        // Check minimum mole fraction
        double minMolFrac = thermo.dMolFraction[iFirst];
        for (int i = iFirst + 1; i < iLast; ++i) {
            minMolFrac = std::min(minMolFrac, thermo.dMolFraction[i]);
        }
        if (minMolFrac < dMinMoleFraction) break;

        // Compute chemical potentials
        subMinChemicalPotential(ctx, iSolnPhaseIndex, iFirst, iLast);

        // Compute driving force
        dDrivingForce = subMinDrivingForce(ctx, iFirst, iLast, dChemicalPotentialStar);

        if (dDrivingForce < dDrivingForceLast) {
            break;  // Wolfe condition satisfied
        } else {
            // Dampen the step
            dStepLength *= 0.5;
            dMaxUpdate = 0.0;

            for (int j = 0; j < nVar; ++j) {
                int i = iFirst + j;
                thermo.dMolFraction[i] -= dStepLength * dRHS(j);
                thermo.dMolFraction[i] = std::max(thermo.dMolFraction[i], dMinMoleFraction);
                dMaxUpdate = std::max(dMaxUpdate, std::abs(dStepLength * dRHS(j)));
            }
        }
    }

    // Check convergence
    if (dMaxUpdate <= dSubMinTolerance) {
        lSubMinConverged = true;
    }

    dDrivingForceLast = dDrivingForce;
}

/// @brief Check if two phases in a miscibility gap have duplicate compositions
static bool subMinCheckDuplicate(ThermoContext& ctx, int /*iSolnPhaseIndex*/,
                                 int iSolnPhaseIndexOther,
                                 int iFirst, int /*iLast*/, int nVar) {
    auto& thermo = *ctx.thermo;

    const double dTolEuclideanNorm = 1e-4;

    int iFirstOther = thermo.nSpeciesPhase[iSolnPhaseIndexOther - 1] + 1;
    int iLastOther = thermo.nSpeciesPhase[iSolnPhaseIndexOther];

    if (iLastOther - iFirstOther + 1 != nVar) {
        return false;  // Different number of constituents
    }

    // Compute Euclidean norm between mole fraction vectors
    double dTemp = 0.0;
    for (int k = 0; k < nVar; ++k) {
        int i = iFirst + k;
        int j = iFirstOther + k;
        double diff = thermo.dMolFraction[i] - thermo.dMolFraction[j];
        dTemp += diff * diff;
    }

    dTemp = std::sqrt(dTemp) / static_cast<double>(nVar);

    return dTemp < dTolEuclideanNorm;
}

/// @brief Perform subminimization for a solution phase
/// @details Determines whether a non-ideal solution phase should be added to
/// the system by finding the composition that minimizes the driving force.
/// @param ctx The thermochimica context
/// @param iSolnPhaseIndex Solution phase index
/// @return True if phase should be added to the system
bool subminimization(ThermoContext& ctx, int iSolnPhaseIndex) {
    auto& thermo = *ctx.thermo;
    auto& gem = *ctx.gem;
    auto& io = *ctx.io;

    const double dMinMoleFraction = 1e-100;
    const int iterSubMax = 100;

    bool lPhasePass = false;
    bool lDuplicate = false;
    bool lSubMinConverged = false;

    // Initialize indices (0-based, half-open interval [iFirst, iLast))
    // nSpeciesPhase is cumulative: for phase p, species are at indices
    // nSpeciesPhase[p-1] to nSpeciesPhase[p]-1 (0-based)
    int iFirst = thermo.nSpeciesPhase[iSolnPhaseIndex - 1];
    int iLast = thermo.nSpeciesPhase[iSolnPhaseIndex];
    int nVar = iLast - iFirst;

    // Initialize chemical potential star vector
    std::vector<double> dChemicalPotentialStar(nVar, 0.0);

    // Initialize mole fractions
    for (int k = 0; k < nVar; ++k) {
        int i = iFirst + k;
        thermo.dMolFraction[i] = std::max(thermo.dMolFraction[i], 1e-15);

        // Compute chemical potentials defined by element potentials
        dChemicalPotentialStar[k] = 0.0;
        for (int j = 1; j <= thermo.nElements; ++j) {
            dChemicalPotentialStar[k] +=
                thermo.dElementPotential[j] * thermo.dStoichSpecies(i, j);
        }
        dChemicalPotentialStar[k] /= static_cast<double>(thermo.iParticlesPerMole[i]);
    }

    // Compute initial chemical potentials
    subMinChemicalPotential(ctx, iSolnPhaseIndex, iFirst, iLast);

    // Compute initial driving force
    double dDrivingForce = subMinDrivingForce(ctx, iFirst, iLast, dChemicalPotentialStar);
    double dDrivingForceLast = 10.0;
    double dSubMinFunctionNorm = 0.0;

    // Find other phase for miscibility gap checking
    int iSolnPhaseIndexOther = 0;
    if (gem.lMiscibility[iSolnPhaseIndex]) {
        if (thermo.cSolnPhaseName[iSolnPhaseIndex] == thermo.cSolnPhaseName[iSolnPhaseIndex - 1]) {
            iSolnPhaseIndexOther = iSolnPhaseIndex - 1;
        } else {
            for (int i = 1; i <= thermo.nSolnPhasesSys; ++i) {
                if (i != iSolnPhaseIndex &&
                    thermo.cSolnPhaseName[iSolnPhaseIndex] == thermo.cSolnPhaseName[i]) {
                    iSolnPhaseIndexOther = i;
                    break;
                }
            }
        }
    }

    // Subminimization iteration loop
    Eigen::VectorXd dRHS;
    for (int iterSub = 0; iterSub < iterSubMax; ++iterSub) {
        // Compute Newton direction
        int info = subMinNewton(ctx, iSolnPhaseIndex, iFirst, iLast, nVar,
                               dChemicalPotentialStar, dDrivingForce, dRHS);
        if (info != 0) {
            io.INFOThermo = 28;
            break;
        }

        // Perform line search
        subMinLineSearch(ctx, iSolnPhaseIndex, iFirst, iLast, nVar,
                        dChemicalPotentialStar, dRHS,
                        dDrivingForce, dDrivingForceLast, lSubMinConverged);

        // Compute functional norm
        dSubMinFunctionNorm = subMinFunctionNorm(ctx, iSolnPhaseIndex, iFirst, iLast);

        // Check minimum mole fraction
        double minMolFrac = thermo.dMolFraction[iFirst];
        for (int i = iFirst + 1; i < iLast; ++i) {
            minMolFrac = std::min(minMolFrac, thermo.dMolFraction[i]);
        }
        if (minMolFrac < dMinMoleFraction) break;

        // Check convergence
        if (dSubMinFunctionNorm < thermo.tolerances[kTolMassBalance]) {
            if (dDrivingForce <= thermo.tolerances[kTolDrivingForce] &&
                !gem.lMiscibility[iSolnPhaseIndex]) {
                lSubMinConverged = true;
            }
        }

        if (lSubMinConverged || io.INFOThermo != 0) break;

        // Check for duplicate miscibility gap
        if (gem.lMiscibility[iSolnPhaseIndex] && iSolnPhaseIndexOther > 0) {
            lDuplicate = subMinCheckDuplicate(ctx, iSolnPhaseIndex,
                                              iSolnPhaseIndexOther,
                                              iFirst, iLast, nVar);
            if (lDuplicate) break;
        }
    }

    // Handle results
    if (lDuplicate) {
        dDrivingForce = 0.0;
    }

    if (dSubMinFunctionNorm > 10.0 * thermo.tolerances[kTolMassBalance]) {
        dDrivingForce = 0.0;
    }

    if (dDrivingForce < thermo.tolerances[kTolDrivingForce] && lSubMinConverged) {
        lPhasePass = true;
    }

    gem.dDrivingForceSoln[iSolnPhaseIndex] = dDrivingForce;

    return lPhasePass;
}

/// @brief Compute driving force for a solution phase
double computeDrivingForce(ThermoContext& ctx, int phaseIndex) {
    subminimization(ctx, phaseIndex);
    return ctx.gem->dDrivingForceSoln[phaseIndex];
}

} // namespace Thermochimica
