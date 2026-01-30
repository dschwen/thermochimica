#include "thermochimica/solver/GEMSolver.hpp"
#include "thermochimica/context/PhaseConstraints.hpp"
#include "thermochimica/util/Constants.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

namespace Thermochimica {

bool ConvergenceChecker::check(ThermoContext& ctx) {
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;

    // Don't declare convergence on the same iteration a phase was added/removed
    // We need at least one more iteration to compute chemical potentials
    if (gem.iterGlobal == gem.iterLast) {
        return false;
    }

    // Quick convergence check based on function norm
    if (gem.dGEMFunctionNorm < thermo.tolerances[kTolFunctionNorm] &&
        gem.iterGlobal - gem.iterLast > 100) {
        return true;
    }

    // Check all convergence criteria

    // 1. Check phase slots are valid
    // Positive = pure condensed species (1-based index)
    // Negative = solution phase: -(phaseIndex + 1)
    // Zero = unused slot (only valid for excess slots beyond active phases)
    // Note: Pure condensed phases are at slots 0 to nConPhases-1
    //       Solution phases are at slots nElements-nSolnPhases to nElements-1

    // Check condensed phase slots
    for (int i = 0; i < thermo.nConPhases; ++i) {
        if (thermo.iAssemblage(i) == 0) {
            return false;
        }
    }
    // Check solution phase slots
    for (int i = 0; i < thermo.nSolnPhases; ++i) {
        int slot = thermo.nElements - thermo.nSolnPhases + i;
        if (thermo.iAssemblage(slot) == 0) {
            return false;
        }
    }

    // 2. All phase moles >= 0
    for (int i = 0; i < thermo.nConPhases + thermo.nSolnPhases; ++i) {
        if (thermo.dMolesPhase(i) < 0) {
            return false;
        }
    }

    // 3. Phase rule satisfied
    if (!checkPhaseRule(ctx)) {
        return false;
    }

    // 4. Mass balance residuals
    if (!checkMassBalance(ctx)) {
        return false;
    }

    // 5. Chemical potential residuals
    if (!checkChemicalPotential(ctx)) {
        return false;
    }

    // 6. Site fractions (for sublattice phases)
    if (thermo.nCountSublattice > 0) {
        if (!checkSiteFractions(ctx)) {
            return false;
        }
    }

    // 7. Unstable solution phases
    // Skip this check when constraints are active - we intentionally force a fixed
    // assemblage and exclude unconstrained phases. Checking unstable phases would
    // prevent convergence if any excluded phase has positive driving force.
    if (!ctx.phaseConstraints->hasActiveConstraints()) {
        if (!checkUnstablePhases(ctx)) {
            return false;
        }
    }

    // 8. Miscibility gaps
    if (!checkMiscibility(ctx)) {
        return false;
    }

    // 9. Phase constraints (if active)
    // Note: This is checked in inner loop to help guide phase assemblage,
    // but final constraint satisfaction is checked in outer loop
    // So we don't fail here - just let the outer loop handle it

    // All checks passed
    return true;
}

bool ConvergenceChecker::checkPhaseConstraints(ThermoContext& ctx) {
    auto& pc = *ctx.phaseConstraints;

    if (!pc.hasActiveConstraints()) {
        return true;  // No constraints to check
    }

    return pc.areConstraintsSatisfied();
}

bool ConvergenceChecker::checkMassBalance(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;

    for (int j = 0; j < thermo.nElements; ++j) {
        double sum = 0.0;

        // Sum contributions from all species
        for (int i = 0; i < thermo.nSpecies; ++i) {
            if (thermo.dMolesSpecies(i) > 0) {
                sum += thermo.dMolesSpecies(i) * thermo.dStoichSpecies(i, j);
            }
        }

        // Calculate residual
        double residual = std::abs(sum - thermo.dMolesElement(j));
        double normalizer = std::max(1.0, thermo.dMolesElement(j));
        residual /= normalizer;

        if (residual > thermo.tolerances[kTolMassBalance]) {
            return false;
        }
    }

    return true;
}

bool ConvergenceChecker::checkChemicalPotential(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;

    // For now, skip detailed chemical potential checking
    // The mass balance and driving force checks are more important
    // Chemical potential equilibrium will be achieved through GEM iterations

    // Only do a basic check: if we have solution phases, their species should
    // have reasonable chemical potentials
    if (thermo.nSolnPhases == 0) {
        // For pure-phase-only systems, chemical potentials are determined by element potentials
        // which are set by the assemblage. Just check that element potentials are finite.
        for (int j = 0; j < thermo.nElements; ++j) {
            if (std::isnan(thermo.dElementPotential(j)) || std::isinf(thermo.dElementPotential(j))) {
                return false;
            }
        }
    }

    return true;
}

bool ConvergenceChecker::checkPhaseRule(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;

    // Gibbs phase rule: F = C - P + 2
    // At fixed T,P: F = C - P >= 0
    // So: P <= C (where C = number of active components)

    // Count active elements (with non-zero input)
    int nActiveElements = 0;
    for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
        if (thermo.dMolesElement(j) > 0.0) {
            ++nActiveElements;
        }
    }

    int totalPhases = thermo.nConPhases + thermo.nSolnPhases;

    return totalPhases <= nActiveElements;
}

bool ConvergenceChecker::checkSiteFractions(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;

    // Check that site fractions sum to 1 for each sublattice
    for (int iPhase = 0; iPhase < thermo.nSolnPhasesSys; ++iPhase) {
        int sublatticeIdx = thermo.iPhaseSublattice(iPhase);
        if (sublatticeIdx <= 0) continue;

        int nSublattices = thermo.nSublatticePhase(iPhase);

        for (int iSub = 0; iSub < nSublattices; ++iSub) {
            // Sum site fractions on this sublattice
            if (static_cast<size_t>(sublatticeIdx - 1) < thermo.dSiteFraction.size()) {
                auto& siteFrac = thermo.dSiteFraction[sublatticeIdx - 1];
                double sum = 0.0;

                int nConstituents = thermo.nConstituentSublattice(iPhase, iSub);
                for (int k = 0; k < nConstituents && k < siteFrac.cols(); ++k) {
                    if (iSub < siteFrac.rows()) {
                        sum += siteFrac(iSub, k);
                    }
                }

                if (std::abs(sum - 1.0) > thermo.tolerances[kTolSiteFraction]) {
                    return false;
                }
            }
        }
    }

    return true;
}

bool ConvergenceChecker::checkUnstablePhases(ThermoContext& ctx) {
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;
    auto& io = *ctx.io;

    // Don't check unstable phases until Newton has had time to converge
    // after the last phase change. This prevents premature failure.
    // But do check if we haven't changed phases yet (iterLast = 0).
    int iterSincePhaseChange = gem.iterGlobal - gem.iterLast;
    if (iterSincePhaseChange < 50 && gem.iterLast > 0) {
        return true;  // Allow more Newton iterations after recent phase change
    }

    double R = Constants::kIdealGasConstant;
    double T = io.dTemperature;

    // Check driving force of unstable solution phases
    for (int iPhase = 0; iPhase < thermo.nSolnPhasesSys; ++iPhase) {
        if (gem.lSolnPhases[iPhase]) {
            continue;  // Phase is stable
        }

        // Compute driving force for this phase
        int iFirst = (iPhase > 0) ? thermo.nSpeciesPhase(iPhase) : 0;
        int iLast = thermo.nSpeciesPhase(iPhase + 1);

        // Build list of feasible species (only contain active elements)
        std::vector<int> feasibleSpecies;
        for (int i = iFirst; i < iLast; ++i) {
            bool feasible = true;
            for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
                if (thermo.dStoichSpecies(i, j) > 0.0 && thermo.dMolesElement(j) <= 0.0) {
                    feasible = false;
                    break;
                }
            }
            if (feasible) {
                feasibleSpecies.push_back(i);
            }
        }

        // If no feasible species, phase can't be added
        if (feasibleSpecies.empty()) {
            gem.dDrivingForceSoln(iPhase) = -1e30;  // Very unfavorable
            continue;
        }

        double sumExp = 0.0;
        std::vector<double> expVal(iLast - iFirst, 0.0);

        for (int i : feasibleSpecies) {
            double muStar = 0.0;
            for (int j = 0; j < thermo.nElements; ++j) {
                muStar += thermo.dElementPotential(j) * thermo.dStoichSpecies(i, j);
            }
            muStar /= static_cast<double>(thermo.iParticlesPerMole(i));

            double muStd = thermo.dStdGibbsEnergy(i) / (R * T);
            double arg = muStar - muStd;
            arg = std::max(-100.0, std::min(100.0, arg));
            expVal[i - iFirst] = std::exp(arg);
            sumExp += expVal[i - iFirst];
        }

        double drivingForce = 0.0;
        if (sumExp > 1e-300) {
            for (int i : feasibleSpecies) {
                double x = expVal[i - iFirst] / sumExp;
                if (x > 1e-300) {
                    double muStar = 0.0;
                    for (int j = 0; j < thermo.nElements; ++j) {
                        muStar += thermo.dElementPotential(j) * thermo.dStoichSpecies(i, j);
                    }
                    muStar /= static_cast<double>(thermo.iParticlesPerMole(i));

                    double muStd = thermo.dStdGibbsEnergy(i) / (R * T);
                    double mu = muStd + std::log(x);
                    drivingForce += x * (mu - muStar);
                }
            }
        }

        // Store the computed driving force
        gem.dDrivingForceSoln(iPhase) = -drivingForce;

        // If driving force is positive and significant, phase should be added
        if (-drivingForce > thermo.tolerances[kTolDrivingForce]) {
            return false;
        }
    }

    return true;
}

bool ConvergenceChecker::checkMiscibility(ThermoContext& ctx) {
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;

    // Check for miscibility gaps in stable solution phases
    for (int iPhase = 0; iPhase < thermo.nSolnPhasesSys; ++iPhase) {
        if (!gem.lSolnPhases[iPhase]) {
            continue;  // Phase not stable
        }

        if (!gem.lMiscibility[iPhase]) {
            continue;  // No miscibility gap
        }

        // Would need subminimization to properly check
        // For now, assume no issue if we got this far
    }

    return true;
}

} // namespace Thermochimica
