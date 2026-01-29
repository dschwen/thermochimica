/// @file LevelingSolver.cpp
/// @brief Leveling solver for initial phase assemblage estimation
/// @author Converted from Fortran LevelingSolver.f90

#include "thermochimica/ThermoContext.hpp"
#include "thermochimica/util/Tolerances.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <limits>

namespace Thermochimica {

// Forward declarations
void getFirstAssemblage(ThermoContext& ctx);
void getNewAssemblage(ThermoContext& ctx, int iter);
void computeLevel(ThermoContext& ctx);
static bool speciesIsFeasible(ThermoContext& ctx, int speciesIdx);

/// @brief Linear solver for initial equilibrium estimation
/// @details Uses the "Leveling" technique of Eriksson and Thompson.
/// Temporarily treats all species as pure stoichiometric phases,
/// making the Gibbs energy of the system a linear function.
///
/// Advantages:
/// 1. Provides close estimates for element potentials and chemical potentials
/// 2. Estimates phase assemblage close to equilibrium
/// 3. Determines assemblage with few iterations
/// 4. Provides estimates for mole fractions of solution phase constituents
///
/// References:
/// - G. Eriksson and W.T. Thompson, CALPHAD V.13, N.4, pp.389-400 (1989)
/// - M.H.A. Piro, PhD Dissertation, Royal Military College of Canada (2011)
///
/// @param ctx The thermochimica context
void levelingSolver(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;
    auto& io = *ctx.io;

    int nElements = thermo.nElements;
    int nSpecies = thermo.nSpecies;

    // Allocate/reallocate arrays (using 0-based indexing)
    thermo.dMolesPhase.resize(nElements);
    thermo.dElementPotential.resize(nElements);
    thermo.iAssemblage.resize(nElements);
    thermo.dLevel.resize(nElements);
    thermo.iterHistoryLevel.resize(nElements, 1000);

    // Initialize variables
    thermo.dElementPotential.setZero();
    thermo.iterHistoryLevel.setZero();
    thermo.iAssemblage.setZero();  // 0 means no species assigned
    thermo.dLevel.setZero();
    thermo.dMolesPhase.setZero();
    thermo.nConPhases = nElements;

    // Establish the very first phase assemblage
    getFirstAssemblage(ctx);

    if (io.INFOThermo != 0) return;

    // Build list of active elements (with non-zero input)
    std::vector<int> activeElements;
    for (int e = 0; e < nElements - thermo.nChargedConstraints; ++e) {
        if (thermo.dMolesElement(e) > 0.0) {
            activeElements.push_back(e);
        }
    }

    // Index of first pure condensed species
    int pureStart = thermo.nSpeciesPhase(thermo.nSolnPhasesSys);

    // START LEVELING
    // Limit iterations - if we can't converge with pure phases, the GEM solver
    // will add solution phases as needed
    const int maxIter = 100;
    int iter;

    double R = Constants::kIdealGasConstant;
    double T = io.dTemperature;

    for (iter = 0; iter < maxIter; ++iter) {
        // Compute element potentials for current assemblage
        // (computeLevel now uses dStdGibbsEnergy directly, not modified dChemicalPotential)
        computeLevel(ctx);

        // Calculate driving force for all feasible pure species based on element potentials
        // For species i: DF_i = G°_i/RT - sum_j(lambda_j * a_ij)
        double minDrivingForce = std::numeric_limits<double>::max();

        // Only check feasible pure condensed species
        for (int i = pureStart; i < nSpecies; ++i) {
            // Check if species is feasible (all elements have non-zero input)
            if (!speciesIsFeasible(ctx, i)) continue;

            double muStar = 0.0;  // sum_j(lambda_j * a_ij)
            for (int elemIdx : activeElements) {
                muStar += thermo.dElementPotential(elemIdx) * thermo.dStoichSpecies(i, elemIdx);
            }
            // Driving force = G°/RT - mu* (negative means species is unstable)
            double drivingForce = thermo.dStdGibbsEnergy(i) / (R * T) - muStar;
            minDrivingForce = std::min(minDrivingForce, drivingForce);
        }

        // Leveling is complete when all feasible species have non-negative driving force
        if (minDrivingForce >= -thermo.tolerances[kTolDrivingForce]) {
            break;
        }

        // Determine next phase assemblage
        getNewAssemblage(ctx, iter);

        if (io.INFOThermo != 0) break;
    }

    // Compute moles of phases after leveling
    computeLevel(ctx);

}

/// @brief Check if a species only contains elements that have non-zero input amounts
/// @param ctx The thermochimica context
/// @param speciesIdx The species index (0-based)
/// @return true if species is "feasible" (all its elements have non-zero input)
static bool speciesIsFeasible(ThermoContext& ctx, int speciesIdx) {
    auto& thermo = *ctx.thermo;
    int nElements = thermo.nElements - thermo.nChargedConstraints;

    for (int e = 0; e < nElements; ++e) {
        // If species contains this element but element has zero input, species is infeasible
        if (thermo.dStoichSpecies(speciesIdx, e) > 0.0 &&
            thermo.dMolesElement(e) <= 0.0) {
            return false;
        }
    }
    return true;
}

/// @brief Get first phase assemblage using species with lowest Gibbs energies
void getFirstAssemblage(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;

    int nElements = thermo.nElements;
    int nSpecies = thermo.nSpecies;

    // Index of first pure condensed species
    int pureStart = thermo.nSpeciesPhase(thermo.nSolnPhasesSys);

    // For each element with non-zero input, find feasible species with lowest
    // chemical potential that contains that element (0-based indexing)
    // ONLY consider pure condensed phases (species >= pureStart)

    std::vector<bool> speciesUsed(nSpecies, false);

    // Only assign phases for active elements
    int assignedPhases = 0;
    for (int e = 0; e < nElements - thermo.nChargedConstraints; ++e) {
        // Skip elements with zero input
        if (thermo.dMolesElement(e) <= 0.0) {
            continue;
        }

        double minGibbs = std::numeric_limits<double>::max();
        int bestSpecies = -1;  // -1 means no species assigned

        // Only consider feasible pure condensed phases
        for (int i = pureStart; i < nSpecies; ++i) {
            if (speciesUsed[i]) continue;

            // Check if species is feasible (all elements have non-zero input)
            if (!speciesIsFeasible(ctx, i)) continue;

            // Check if species contains this element
            if (thermo.dStoichSpecies(i, e) > 0.0) {
                double gibbs = thermo.dChemicalPotential(i);
                if (gibbs < minGibbs) {
                    minGibbs = gibbs;
                    bestSpecies = i;
                }
            }
        }

        if (bestSpecies >= 0) {
            thermo.iAssemblage(assignedPhases) = bestSpecies + 1;  // Store as 1-based for compatibility
            speciesUsed[bestSpecies] = true;
            assignedPhases++;
        }
    }

    // Update nConPhases to reflect actual assigned phases
    thermo.nConPhases = assignedPhases;

    // Compute initial level adjustments
    computeLevel(ctx);
}

/// @brief Get new phase assemblage by swapping species
/// @details Find pure condensed species with most negative driving force
/// and swap it with the worst species in the current assemblage
void getNewAssemblage(ThermoContext& ctx, int iter) {
    auto& thermo = *ctx.thermo;
    auto& io = *ctx.io;

    int nElements = thermo.nElements;
    int nSpecies = thermo.nSpecies;
    int nConPhases = thermo.nConPhases;
    double R = Constants::kIdealGasConstant;
    double T = io.dTemperature;

    // Index of first pure condensed species
    int pureStart = thermo.nSpeciesPhase(thermo.nSolnPhasesSys);

    // Find feasible pure condensed species with most negative driving force
    // Driving force = G°/RT - sum_j(lambda_j * a_ij)
    // Negative driving force means species is more stable than current potential allows
    int iMinSpecies = -1;
    double minDrivingForce = 0.0;

    // Only consider feasible pure condensed phases
    for (int i = pureStart; i < nSpecies; ++i) {
        // Check if species is feasible (all elements have non-zero input)
        if (!speciesIsFeasible(ctx, i)) continue;

        // Compute driving force
        double muStar = 0.0;
        for (int j = 0; j < nElements - thermo.nChargedConstraints; ++j) {
            muStar += thermo.dElementPotential(j) * thermo.dStoichSpecies(i, j);
        }
        double drivingForce = thermo.dStdGibbsEnergy(i) / (R * T) - muStar;

        if (drivingForce < minDrivingForce) {
            minDrivingForce = drivingForce;
            iMinSpecies = i;
        }
    }

    if (iMinSpecies < 0) {
        return;  // All feasible pure species have non-negative driving force
    }

    // Check if this species is already in the assemblage (stored as 1-based)
    for (int p = 0; p < nConPhases; ++p) {
        if (thermo.iAssemblage(p) == iMinSpecies + 1) {
            return;  // Already in assemblage
        }
    }

    // Find best position to add this species
    // Replace the species that has the most positive driving force (least stable)
    double bestImprovement = 0.0;
    int bestPosition = -1;

    // Build list of active elements for computing mu*
    std::vector<int> activeElements;
    for (int e = 0; e < nElements - thermo.nChargedConstraints; ++e) {
        if (thermo.dMolesElement(e) > 0.0) {
            activeElements.push_back(e);
        }
    }

    for (int p = 0; p < nConPhases; ++p) {
        int currentSpecies = thermo.iAssemblage(p) - 1;  // Convert to 0-based
        if (currentSpecies < 0) continue;

        // Check if new species shares elements with current species
        bool sharesElement = false;
        for (int elemIdx : activeElements) {
            if (thermo.dStoichSpecies(iMinSpecies, elemIdx) > 0.0 &&
                thermo.dStoichSpecies(currentSpecies, elemIdx) > 0.0) {
                sharesElement = true;
                break;
            }
        }
        if (!sharesElement) continue;

        // Compute driving force for current species
        double muStarCurrent = 0.0;
        for (int elemIdx : activeElements) {
            muStarCurrent += thermo.dElementPotential(elemIdx) * thermo.dStoichSpecies(currentSpecies, elemIdx);
        }
        double dfCurrent = thermo.dStdGibbsEnergy(currentSpecies) / (R * T) - muStarCurrent;

        // Improvement = how much better the new species is
        double improvement = dfCurrent - minDrivingForce;
        if (improvement > bestImprovement) {
            bestImprovement = improvement;
            bestPosition = p;
        }
    }

    // Perform the swap
    if (bestPosition >= 0) {
        // Record history to prevent cycling
        thermo.iterHistoryLevel(bestPosition, iter % 1000) = iMinSpecies + 1;

        // Check for cycling
        bool cycling = false;
        for (int prev = 0; prev < iter && prev < 1000; ++prev) {
            if (thermo.iterHistoryLevel(bestPosition, prev) == iMinSpecies + 1) {
                cycling = true;
                break;
            }
        }

        if (!cycling) {
            thermo.iAssemblage(bestPosition) = iMinSpecies + 1;
        }
    }
}

/// @brief Compute level adjustments and element potentials based on current assemblage
/// @details Computes element potentials that would give zero chemical potential
/// for the species in the current assemblage. For equilibrium:
/// mu_i = G°_i/RT = sum_j(lambda_j * a_ij)  =>  lambda = A^-T * (G°/RT)
void computeLevel(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;
    auto& io = *ctx.io;

    int nElements = thermo.nElements;
    int nConPhases = thermo.nConPhases;  // Use actual assigned phases, not nElements
    double R = Constants::kIdealGasConstant;
    double T = io.dTemperature;

    if (nConPhases == 0) return;

    // Build list of active elements (with non-zero input)
    std::vector<int> activeElements;
    for (int e = 0; e < nElements - thermo.nChargedConstraints; ++e) {
        if (thermo.dMolesElement(e) > 0.0) {
            activeElements.push_back(e);
        }
    }
    int nActiveElements = activeElements.size();

    // If nConPhases != nActiveElements, we have a mismatch
    // Use min to build a solvable system
    int systemSize = std::min(nConPhases, nActiveElements);
    if (systemSize == 0) return;

    // Build stoichiometry matrix from current assemblage (0-based indexing)
    // Only for active elements
    Eigen::MatrixXd stoichMatrix = Eigen::MatrixXd::Zero(systemSize, systemSize);
    Eigen::VectorXd gibbsVec = Eigen::VectorXd::Zero(systemSize);

    for (int p = 0; p < systemSize; ++p) {
        int species = thermo.iAssemblage(p) - 1;  // Convert from 1-based to 0-based
        if (species < 0) continue;

        // Use standard Gibbs energy normalized by RT (not the modified dChemicalPotential)
        gibbsVec(p) = thermo.dStdGibbsEnergy(species) / (R * T);

        for (int j = 0; j < systemSize; ++j) {
            int elemIdx = activeElements[j];
            stoichMatrix(p, j) = thermo.dStoichSpecies(species, elemIdx);
        }
    }

    // Check if matrix is well-conditioned
    double det = stoichMatrix.determinant();

    if (std::abs(det) < 1e-20) {
        return;  // Singular or near-singular matrix
    }

    // Solve for element potentials directly
    // For pure phase at equilibrium: G°/RT = sum_j(lambda_j * a_ij)
    // In matrix form: A * lambda = g°  =>  lambda = A^-1 * g°
    // stoichMatrix rows correspond to phases, cols to elements
    // stoichMatrix * lambda = gibbs => lambda = stoichMatrix^-1 * gibbs
    Eigen::VectorXd lambda = stoichMatrix.colPivHouseholderQr().solve(gibbsVec);

    // Set element potentials for active elements only
    // Map lambda back to original element indices
    for (int j = 0; j < systemSize; ++j) {
        int elemIdx = activeElements[j];
        thermo.dElementPotential(elemIdx) = lambda(j);
        thermo.dLevel(elemIdx) = lambda(j);  // Store for reference
    }

    // Compute moles of phases from mass balance:
    // sum_p (n_p * a_pj) = b_j  for each element j
    // In matrix form: A^T * n = b where A is (phases x elements)
    // So: stoichMatrix^T * molesPhase = molesElement
    Eigen::VectorXd molesElement(systemSize);
    for (int j = 0; j < systemSize; ++j) {
        int elemIdx = activeElements[j];
        molesElement(j) = thermo.dMolesElement(elemIdx);
    }

    Eigen::VectorXd molesPhase = stoichMatrix.transpose().colPivHouseholderQr().solve(molesElement);

    for (int p = 0; p < systemSize; ++p) {
        thermo.dMolesPhase(p) = std::max(molesPhase(p), 0.0);

        // Also set dMolesSpecies for pure condensed phases
        int species = thermo.iAssemblage(p) - 1;  // Convert from 1-based to 0-based
        if (species >= 0 && species < thermo.nSpecies) {
            thermo.dMolesSpecies(species) = thermo.dMolesPhase(p);
            thermo.dMolFraction(species) = 1.0;  // Pure phase has mole fraction = 1
        }
    }
}

/// @brief Post-leveling adjustments
void postLevelingSolver(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;
    auto& io = *ctx.io;
    double R = Constants::kIdealGasConstant;
    double T = io.dTemperature;

    // Compute initial mole fractions for solution phases (0-based indexing)
    for (int p = 0; p < thermo.nSolnPhasesSys; ++p) {
        int iFirst = (p > 0) ? thermo.nSpeciesPhase(p) : 0;
        int iLast = thermo.nSpeciesPhase(p + 1);

        // Estimate mole fractions from element potentials
        // x_i ∝ exp(mu* - mu_std) where both are dimensionless
        double sumExp = 0.0;
        for (int i = iFirst; i < iLast; ++i) {
            // mu* = sum_j(lambda_j * a_ij) / n_i  (dimensionless)
            double muStar = 0.0;
            for (int j = 0; j < thermo.nElements; ++j) {
                muStar += thermo.dElementPotential(j) * thermo.dStoichSpecies(i, j);
            }
            muStar /= static_cast<double>(thermo.iParticlesPerMole(i));

            // mu_std = G°/RT (dimensionless)
            double muStd = thermo.dStdGibbsEnergy(i) / (R * T);
            double arg = muStar - muStd;
            arg = std::max(-100.0, std::min(100.0, arg));  // Prevent overflow
            double expVal = std::exp(arg);
            thermo.dMolFraction(i) = expVal;
            sumExp += expVal;
        }

        // Normalize
        if (sumExp > 0.0) {
            for (int i = iFirst; i < iLast; ++i) {
                thermo.dMolFraction(i) /= sumExp;
                thermo.dMolFraction(i) = std::max(thermo.dMolFraction(i), 1e-100);
            }
        }
    }
}

} // namespace Thermochimica
