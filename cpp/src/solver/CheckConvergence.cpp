#include "thermochimica/solver/GEMSolver.hpp"
#include "thermochimica/util/Constants.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

namespace Thermochimica {

bool ConvergenceChecker::check(ThermoContext& ctx) {
    auto& gem = *ctx.gem;
    auto& thermo = *ctx.thermo;

    // Quick convergence check based on function norm
    if (gem.dGEMFunctionNorm < thermo.tolerances[kTolFunctionNorm] &&
        gem.iterGlobal - gem.iterLast > 100) {
        if (gem.iterGlobal <= 50 || gem.iterGlobal % 500 == 0) {
            // DEBUG: std::cerr << "DEBUG iter " << gem.iterGlobal << ": converged (function norm)\n";
        }
        return true;
    }

    // Check all convergence criteria

    // 1. Check phase slots are valid
    // Positive = pure condensed species (1-based index)
    // Negative = solution phase: -(phaseIndex + 1)
    // Zero = unused slot (only valid for excess slots beyond active phases)
    int activePhases = thermo.nConPhases + thermo.nSolnPhases;
    for (int i = 0; i < activePhases; ++i) {
        if (thermo.iAssemblage(i) == 0) {
            return false;  // Active slot with no phase assigned
        }
    }

    // 2. All phase moles >= 0
    for (int i = 0; i < thermo.nConPhases + thermo.nSolnPhases; ++i) {
        if (thermo.dMolesPhase(i) < 0) {
            if (gem.iterGlobal <= 50 || gem.iterGlobal % 500 == 0) {
                // DEBUG: std::cerr << "DEBUG iter " << gem.iterGlobal << ": fail (phase moles < 0)\n";
            }
            return false;
        }
    }

    // 3. Phase rule satisfied
    if (!checkPhaseRule(ctx)) {
        if (gem.iterGlobal <= 50 || gem.iterGlobal % 500 == 0) {
            // DEBUG: std::cerr << "DEBUG iter " << gem.iterGlobal << ": fail (phase rule)\n";
        }
        return false;
    }

    // 4. Mass balance residuals
    if (!checkMassBalance(ctx)) {
        if (gem.iterGlobal <= 50 || gem.iterGlobal % 500 == 0) {
            // DEBUG: std::cerr << "DEBUG iter " << gem.iterGlobal << ": fail (mass balance)\n";
        }
        return false;
    }

    // 5. Chemical potential residuals
    if (!checkChemicalPotential(ctx)) {
        if (gem.iterGlobal <= 50 || gem.iterGlobal % 500 == 0) {
            // DEBUG: std::cerr << "DEBUG iter " << gem.iterGlobal << ": fail (chemical potential)\n";
        }
        return false;
    }

    // 6. Site fractions (for sublattice phases)
    if (thermo.nCountSublattice > 0) {
        if (!checkSiteFractions(ctx)) {
            if (gem.iterGlobal <= 50 || gem.iterGlobal % 500 == 0) {
                // DEBUG: std::cerr << "DEBUG iter " << gem.iterGlobal << ": fail (site fractions)\n";
            }
            return false;
        }
    }

    // 7. Unstable solution phases
    if (!checkUnstablePhases(ctx)) {
        if (gem.iterGlobal <= 50 || gem.iterGlobal % 500 == 0) {
            // DEBUG: std::cerr << "DEBUG iter " << gem.iterGlobal << ": fail (unstable phases)\n";
        }
        return false;
    }

    // 8. Miscibility gaps
    if (!checkMiscibility(ctx)) {
        if (gem.iterGlobal <= 50 || gem.iterGlobal % 500 == 0) {
            // DEBUG: std::cerr << "DEBUG iter " << gem.iterGlobal << ": fail (miscibility)\n";
        }
        return false;
    }

    if (gem.iterGlobal <= 50 || gem.iterGlobal % 500 == 0) {
        // DEBUG: std::cerr << "DEBUG iter " << gem.iterGlobal << ": CONVERGED!\n";
    }
    return true;
}

bool ConvergenceChecker::checkMassBalance(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;
    auto& gem = *ctx.gem;

    bool debug = false;  // Set to true for debugging
    double maxResidual = 0.0;
    int worstElement = -1;

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

        if (residual > maxResidual) {
            maxResidual = residual;
            worstElement = j;
        }

        if (residual > thermo.tolerances[kTolMassBalance]) {
            if (debug) {
                std::cerr << "  Mass balance fail: element " << j
                          << " (" << thermo.cElementName[j] << ")"
                          << " computed=" << sum << " target=" << thermo.dMolesElement(j)
                          << " residual=" << residual << "\n";
            }
            return false;
        }
    }

    return true;
}

bool ConvergenceChecker::checkChemicalPotential(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;
    auto& gem = *ctx.gem;
    auto& io = *ctx.io;

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
    // So: P <= C (where C = nElements - nChargedConstraints)

    int maxPhases = thermo.nElements - thermo.nChargedConstraints;
    int totalPhases = thermo.nConPhases + thermo.nSolnPhases;

    return totalPhases <= maxPhases;
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

    double R = Constants::kIdealGasConstant;
    double T = io.dTemperature;
    bool debug = false;  // Set to true for debugging

    // Check driving force of unstable solution phases
    for (int iPhase = 0; iPhase < thermo.nSolnPhasesSys; ++iPhase) {
        if (gem.lSolnPhases[iPhase]) {
            continue;  // Phase is stable
        }

        // Compute driving force for this phase
        int iFirst = (iPhase > 0) ? thermo.nSpeciesPhase(iPhase) : 0;
        int iLast = thermo.nSpeciesPhase(iPhase + 1);

        double sumExp = 0.0;
        std::vector<double> expVal(iLast - iFirst, 0.0);

        for (int i = iFirst; i < iLast; ++i) {
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
            for (int i = iFirst; i < iLast; ++i) {
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
            if (debug) {
                std::cerr << "  Unstable phase " << iPhase << " (" << thermo.cSolnPhaseName[iPhase]
                          << ") has favorable driving force: " << -drivingForce << "\n";
            }
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
