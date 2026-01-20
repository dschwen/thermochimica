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

    // START LEVELING
    // Limit iterations - if we can't converge with pure phases, the GEM solver
    // will add solution phases as needed
    const int maxIter = 100;
    int iter;

    for (iter = 0; iter < maxIter; ++iter) {
        // Compute element potentials for current assemblage
        // (computeLevel now uses dStdGibbsEnergy directly, not modified dChemicalPotential)
        computeLevel(ctx);

        // Calculate chemical potentials for all species based on element potentials
        // For species i: mu_i = G°_i/RT - sum_j(lambda_j * a_ij)
        double minChemPot = std::numeric_limits<double>::max();
        double R = Constants::kIdealGasConstant;
        double T = io.dTemperature;

        for (int i = 0; i < nSpecies; ++i) {
            double muStar = 0.0;  // sum_j(lambda_j * a_ij)
            for (int j = 0; j < nElements; ++j) {
                muStar += thermo.dElementPotential(j) * thermo.dStoichSpecies(i, j);
            }
            // Driving force = G°/RT - mu* (negative means species is unstable)
            double drivingForce = thermo.dStdGibbsEnergy(i) / (R * T) - muStar;
            minChemPot = std::min(minChemPot, drivingForce);
        }

        // Leveling is complete when all species have non-negative driving force
        if (minChemPot >= -thermo.tolerances[kTolDrivingForce]) {
            break;
        }

        // Determine next phase assemblage
        getNewAssemblage(ctx, iter);

        if (io.INFOThermo != 0) break;
    }

    // Compute moles of phases after leveling
    computeLevel(ctx);

}

/// @brief Get first phase assemblage using species with lowest Gibbs energies
void getFirstAssemblage(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;

    int nElements = thermo.nElements;
    int nSpecies = thermo.nSpecies;

    // Index of first pure condensed species
    int pureStart = thermo.nSpeciesPhase(thermo.nSolnPhasesSys);

    // For each element, find species with lowest chemical potential
    // that contains that element (0-based indexing)
    // ONLY consider pure condensed phases (species >= pureStart)

    std::vector<bool> speciesUsed(nSpecies, false);

    for (int e = 0; e < nElements - thermo.nChargedConstraints; ++e) {
        double minGibbs = std::numeric_limits<double>::max();
        int bestSpecies = -1;  // -1 means no species assigned

        // Only consider pure condensed phases
        for (int i = pureStart; i < nSpecies; ++i) {
            if (speciesUsed[i]) continue;

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
            thermo.iAssemblage(e) = bestSpecies + 1;  // Store as 1-based for compatibility
            speciesUsed[bestSpecies] = true;
        }
    }

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
    double R = Constants::kIdealGasConstant;
    double T = io.dTemperature;

    // Index of first pure condensed species
    int pureStart = thermo.nSpeciesPhase(thermo.nSolnPhasesSys);

    // Find pure condensed species with most negative driving force
    // Driving force = G°/RT - sum_j(lambda_j * a_ij)
    // Negative driving force means species is more stable than current potential allows
    int iMinSpecies = -1;
    double minDrivingForce = 0.0;

    // Only consider pure condensed phases
    for (int i = pureStart; i < nSpecies; ++i) {
        // Compute driving force
        double muStar = 0.0;
        for (int j = 0; j < nElements; ++j) {
            muStar += thermo.dElementPotential(j) * thermo.dStoichSpecies(i, j);
        }
        double drivingForce = thermo.dStdGibbsEnergy(i) / (R * T) - muStar;

        if (drivingForce < minDrivingForce) {
            minDrivingForce = drivingForce;
            iMinSpecies = i;
        }
    }

    if (iMinSpecies < 0) {
        return;  // All pure species have non-negative driving force
    }

    // Check if this species is already in the assemblage (stored as 1-based)
    for (int e = 0; e < nElements - thermo.nChargedConstraints; ++e) {
        if (thermo.iAssemblage(e) == iMinSpecies + 1) {
            return;  // Already in assemblage
        }
    }

    // Find best position to add this species
    // Replace the species that has the most positive driving force (least stable)
    double bestImprovement = 0.0;
    int bestPosition = -1;

    for (int e = 0; e < nElements - thermo.nChargedConstraints; ++e) {
        int currentSpecies = thermo.iAssemblage(e) - 1;  // Convert to 0-based
        if (currentSpecies < 0) continue;

        // Check if new species covers this element
        if (thermo.dStoichSpecies(iMinSpecies, e) > 0.0) {
            // Compute driving force for current species
            double muStarCurrent = 0.0;
            for (int j = 0; j < nElements; ++j) {
                muStarCurrent += thermo.dElementPotential(j) * thermo.dStoichSpecies(currentSpecies, j);
            }
            double dfCurrent = thermo.dStdGibbsEnergy(currentSpecies) / (R * T) - muStarCurrent;

            // Improvement = how much better the new species is
            double improvement = dfCurrent - minDrivingForce;
            if (improvement > bestImprovement) {
                bestImprovement = improvement;
                bestPosition = e;
            }
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
    int nConPhases = nElements - thermo.nChargedConstraints;
    double R = Constants::kIdealGasConstant;
    double T = io.dTemperature;

    // Build stoichiometry matrix from current assemblage (0-based indexing)
    Eigen::MatrixXd stoichMatrix = Eigen::MatrixXd::Zero(nConPhases, nConPhases);
    Eigen::VectorXd gibbsVec = Eigen::VectorXd::Zero(nConPhases);

    for (int e = 0; e < nConPhases; ++e) {
        int species = thermo.iAssemblage(e) - 1;  // Convert from 1-based to 0-based
        if (species < 0) continue;

        // Use standard Gibbs energy normalized by RT (not the modified dChemicalPotential)
        gibbsVec(e) = thermo.dStdGibbsEnergy(species) / (R * T);

        for (int j = 0; j < nConPhases; ++j) {
            stoichMatrix(e, j) = thermo.dStoichSpecies(species, j);
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
    // But stoichMatrix rows correspond to phases, cols to elements
    // So we need: stoichMatrix^T * lambda = gibbs => lambda = (A^T)^-1 * gibbs
    Eigen::VectorXd lambda = stoichMatrix.transpose().colPivHouseholderQr().solve(gibbsVec);

    // Set element potentials directly (not accumulate)
    for (int j = 0; j < nConPhases; ++j) {
        thermo.dElementPotential(j) = lambda(j);
        thermo.dLevel(j) = lambda(j);  // Store for reference
    }

    // Compute moles of phases
    Eigen::VectorXd molesElement(nConPhases);
    for (int e = 0; e < nConPhases; ++e) {
        molesElement(e) = thermo.dMolesElement(e);
    }

    Eigen::VectorXd molesPhase = stoichMatrix.colPivHouseholderQr().solve(molesElement);

    for (int e = 0; e < nConPhases; ++e) {
        thermo.dMolesPhase(e) = std::max(molesPhase(e), 0.0);

        // Also set dMolesSpecies for pure condensed phases
        int species = thermo.iAssemblage(e) - 1;  // Convert from 1-based to 0-based
        if (species >= 0 && species < thermo.nSpecies) {
            thermo.dMolesSpecies(species) = thermo.dMolesPhase(e);
            thermo.dMolFraction(species) = 1.0;  // Pure phase has mole fraction = 1
        }
    }
}

/// @brief Post-leveling adjustments
void postLevelingSolver(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;

    // Compute initial mole fractions for solution phases (0-based indexing)
    for (int p = 0; p < thermo.nSolnPhasesSys; ++p) {
        int iFirst = (p > 0) ? thermo.nSpeciesPhase(p) : 0;
        int iLast = thermo.nSpeciesPhase(p + 1);

        // Estimate mole fractions from element potentials
        double sumExp = 0.0;
        for (int i = iFirst; i < iLast; ++i) {
            double muStar = 0.0;
            for (int j = 0; j < thermo.nElements; ++j) {
                muStar += thermo.dElementPotential(j) * thermo.dStoichSpecies(i, j);
            }
            muStar /= static_cast<double>(thermo.iParticlesPerMole(i));

            double chemPot = thermo.dStdGibbsEnergy(i);
            double expVal = std::exp(std::min(muStar - chemPot, 100.0));
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
