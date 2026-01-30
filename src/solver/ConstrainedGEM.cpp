#include "thermochimica/solver/GEMSolver.hpp"
#include "thermochimica/context/PhaseConstraints.hpp"
#include "thermochimica/util/Constants.hpp"
#include "thermochimica/util/Tolerances.hpp"
#include <cmath>
#include <algorithm>

namespace Thermochimica {

void ConstrainedGEM::computePhaseElementFractions(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;
    auto& pc = *ctx.phaseConstraints;

    // Compute total element moles in system
    double totalElementMoles = 0.0;
    for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
        totalElementMoles += thermo.dMolesElement(j);
    }

    if (totalElementMoles <= 0.0) {
        // Set all fractions to 0
        for (auto& c : pc.solnPhaseConstraints) {
            c.currentFraction = 0.0;
        }
        for (auto& c : pc.condPhaseConstraints) {
            c.currentFraction = 0.0;
        }
        return;
    }

    // Compute element fraction for each solution phase
    for (int iPhase = 0; iPhase < thermo.nSolnPhasesSys; ++iPhase) {
        if (iPhase >= static_cast<int>(pc.solnPhaseConstraints.size())) {
            continue;
        }

        double phaseElementMoles = 0.0;
        int iFirst = (iPhase > 0) ? thermo.nSpeciesPhase(iPhase) : 0;
        int iLast = thermo.nSpeciesPhase(iPhase + 1);

        for (int i = iFirst; i < iLast; ++i) {
            for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
                phaseElementMoles += thermo.dMolesSpecies(i) * thermo.dStoichSpecies(i, j);
            }
        }

        pc.solnPhaseConstraints[iPhase].currentFraction = phaseElementMoles / totalElementMoles;
    }

    // Compute element fraction for each pure condensed species
    int nSolnSpecies = thermo.nSpeciesPhase(thermo.nSolnPhasesSys);
    for (int iCond = 0; iCond < thermo.nConPhasesSys; ++iCond) {
        if (iCond >= static_cast<int>(pc.condPhaseConstraints.size())) {
            continue;
        }

        int speciesIdx = nSolnSpecies + iCond;
        if (speciesIdx >= thermo.nSpecies) {
            continue;
        }

        double speciesElementMoles = 0.0;
        for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
            speciesElementMoles += thermo.dMolesSpecies(speciesIdx) * thermo.dStoichSpecies(speciesIdx, j);
        }

        pc.condPhaseConstraints[iCond].currentFraction = speciesElementMoles / totalElementMoles;
    }
}

void ConstrainedGEM::updateLagrangeMultipliers(ThermoContext& ctx) {
    auto& pc = *ctx.phaseConstraints;

    // Update Lagrange multipliers: λ += ρ * (f - f_target)
    for (auto& c : pc.solnPhaseConstraints) {
        if (c.mode == PhaseConstraintMode::Fixed) {
            double violation = c.currentFraction - c.targetFraction;
            c.lagrangeMultiplier += pc.penaltyParameter * violation;
        }
    }

    for (auto& c : pc.condPhaseConstraints) {
        if (c.mode == PhaseConstraintMode::Fixed) {
            double violation = c.currentFraction - c.targetFraction;
            c.lagrangeMultiplier += pc.penaltyParameter * violation;
        }
    }
}

void ConstrainedGEM::addConstraintGradients(ThermoContext& ctx,
                                            Eigen::VectorXd& rhs,
                                            const std::vector<int>& activeElements) {
    auto& thermo = *ctx.thermo;
    auto& pc = *ctx.phaseConstraints;

    // Compute total element moles
    double totalElementMoles = 0.0;
    for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
        totalElementMoles += thermo.dMolesElement(j);
    }
    if (totalElementMoles <= 0.0) return;

    int nActiveEl = static_cast<int>(activeElements.size());

    // For each constrained solution phase, add gradient contribution
    // The Lagrangian is: L = G + Σₚ [λₚ(fₚ - fₚ*) + (ρ/2)(fₚ - fₚ*)²]
    // Gradient w.r.t. species moles nᵢ:
    //   ∂L/∂nᵢ = μᵢ + Σₚ [λₚ + ρ(fₚ - fₚ*)] * ∂fₚ/∂nᵢ
    // where ∂fₚ/∂nᵢ = Σⱼ aᵢⱼ / totalElementMoles for species i in phase p

    for (int iPhase = 0; iPhase < thermo.nSolnPhasesSys; ++iPhase) {
        if (iPhase >= static_cast<int>(pc.solnPhaseConstraints.size())) continue;
        auto& c = pc.solnPhaseConstraints[iPhase];
        if (c.mode != PhaseConstraintMode::Fixed) continue;

        double violation = c.currentFraction - c.targetFraction;
        double penaltyFactor = c.lagrangeMultiplier + pc.penaltyParameter * violation;

        int iFirst = (iPhase > 0) ? thermo.nSpeciesPhase(iPhase) : 0;
        int iLast = thermo.nSpeciesPhase(iPhase + 1);

        // Add gradient contribution to element equations (mass balance residuals)
        // The penalty modifies the effective Gibbs energy, which changes element potentials
        for (int jRed = 0; jRed < nActiveEl; ++jRed) {
            int jFull = activeElements[jRed];

            double gradContrib = 0.0;
            for (int i = iFirst; i < iLast; ++i) {
                if (thermo.dMolFraction(i) > 0) {
                    // Contribution from species i to element j
                    double stoichSum = 0.0;
                    for (int k = 0; k < thermo.nElements - thermo.nChargedConstraints; ++k) {
                        stoichSum += thermo.dStoichSpecies(i, k);
                    }
                    gradContrib += thermo.dMolFraction(i) * thermo.dStoichSpecies(i, jFull) *
                                   penaltyFactor / totalElementMoles;
                }
            }
            rhs(jRed) -= gradContrib;
        }
    }

    // For constrained pure condensed phases
    int nSolnSpecies = thermo.nSpeciesPhase(thermo.nSolnPhasesSys);
    for (int iCond = 0; iCond < thermo.nConPhasesSys; ++iCond) {
        if (iCond >= static_cast<int>(pc.condPhaseConstraints.size())) continue;
        auto& c = pc.condPhaseConstraints[iCond];
        if (c.mode != PhaseConstraintMode::Fixed) continue;

        double violation = c.currentFraction - c.targetFraction;
        double penaltyFactor = c.lagrangeMultiplier + pc.penaltyParameter * violation;

        int speciesIdx = nSolnSpecies + iCond;
        if (speciesIdx >= thermo.nSpecies) continue;

        // Add gradient contribution
        for (int jRed = 0; jRed < nActiveEl; ++jRed) {
            int jFull = activeElements[jRed];
            double gradContrib = thermo.dStoichSpecies(speciesIdx, jFull) *
                                 penaltyFactor / totalElementMoles;
            rhs(jRed) -= gradContrib * thermo.dMolesSpecies(speciesIdx);
        }
    }
}

void ConstrainedGEM::addConstraintHessian(ThermoContext& ctx,
                                          Eigen::MatrixXd& hessian,
                                          const std::vector<int>& activeElements) {
    auto& thermo = *ctx.thermo;
    auto& pc = *ctx.phaseConstraints;

    // Compute total element moles
    double totalElementMoles = 0.0;
    for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
        totalElementMoles += thermo.dMolesElement(j);
    }
    if (totalElementMoles <= 0.0) return;

    double invTotal = 1.0 / totalElementMoles;
    double invTotalSq = invTotal * invTotal;
    int nActiveEl = static_cast<int>(activeElements.size());

    // For each constrained solution phase, add Hessian contribution
    // The penalty term adds: (ρ/2)(fₚ - fₚ*)²
    // Hessian contribution: ρ * ∂fₚ/∂nᵢ * ∂fₚ/∂nⱼ

    for (int iPhase = 0; iPhase < thermo.nSolnPhasesSys; ++iPhase) {
        if (iPhase >= static_cast<int>(pc.solnPhaseConstraints.size())) continue;
        auto& c = pc.solnPhaseConstraints[iPhase];
        if (c.mode != PhaseConstraintMode::Fixed) continue;

        int iFirst = (iPhase > 0) ? thermo.nSpeciesPhase(iPhase) : 0;
        int iLast = thermo.nSpeciesPhase(iPhase + 1);

        // Compute effective stoichiometry vector for this phase
        // ∂fₚ/∂λⱼ relates to how element potential changes affect phase fraction
        Eigen::VectorXd dfdLambda(nActiveEl);
        dfdLambda.setZero();

        for (int i = iFirst; i < iLast; ++i) {
            if (thermo.dMolesSpecies(i) <= 0) continue;

            for (int jRed = 0; jRed < nActiveEl; ++jRed) {
                int jFull = activeElements[jRed];
                // Derivative of phase element moles w.r.t. species moles
                double stoichSum = 0.0;
                for (int k = 0; k < thermo.nElements - thermo.nChargedConstraints; ++k) {
                    stoichSum += thermo.dStoichSpecies(i, k);
                }
                dfdLambda(jRed) += thermo.dStoichSpecies(i, jFull) * stoichSum * invTotal;
            }
        }

        // Add outer product contribution to Hessian: ρ * dfdLambda * dfdLambda^T
        for (int iRed = 0; iRed < nActiveEl; ++iRed) {
            for (int jRed = iRed; jRed < nActiveEl; ++jRed) {
                double contrib = pc.penaltyParameter * dfdLambda(iRed) * dfdLambda(jRed) * invTotalSq;
                hessian(iRed, jRed) += contrib;
                if (iRed != jRed) {
                    hessian(jRed, iRed) += contrib;  // Symmetry
                }
            }
        }
    }

    // For constrained pure condensed phases
    int nSolnSpecies = thermo.nSpeciesPhase(thermo.nSolnPhasesSys);
    for (int iCond = 0; iCond < thermo.nConPhasesSys; ++iCond) {
        if (iCond >= static_cast<int>(pc.condPhaseConstraints.size())) continue;
        auto& c = pc.condPhaseConstraints[iCond];
        if (c.mode != PhaseConstraintMode::Fixed) continue;

        int speciesIdx = nSolnSpecies + iCond;
        if (speciesIdx >= thermo.nSpecies) continue;

        // Compute gradient of phase fraction w.r.t. element potentials
        Eigen::VectorXd dfdLambda(nActiveEl);
        dfdLambda.setZero();

        double stoichSum = 0.0;
        for (int k = 0; k < thermo.nElements - thermo.nChargedConstraints; ++k) {
            stoichSum += thermo.dStoichSpecies(speciesIdx, k);
        }

        for (int jRed = 0; jRed < nActiveEl; ++jRed) {
            int jFull = activeElements[jRed];
            dfdLambda(jRed) = thermo.dStoichSpecies(speciesIdx, jFull) * stoichSum * invTotal;
        }

        // Add outer product contribution
        for (int iRed = 0; iRed < nActiveEl; ++iRed) {
            for (int jRed = iRed; jRed < nActiveEl; ++jRed) {
                double contrib = pc.penaltyParameter * dfdLambda(iRed) * dfdLambda(jRed) * invTotalSq;
                hessian(iRed, jRed) += contrib;
                if (iRed != jRed) {
                    hessian(jRed, iRed) += contrib;
                }
            }
        }
    }
}

void ConstrainedGEM::runConstrainedGEMIteration(ThermoContext& ctx) {
    // This function runs a standard GEM iteration with constraint penalty terms
    // The actual iteration logic is in GEMSolver; this function is called
    // to handle constraint-specific updates within the iteration

    // Compute current phase fractions
    computePhaseElementFractions(ctx);
}

bool ConstrainedGEM::setupAssemblageFromConstraints(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;
    auto& gem = *ctx.gem;
    auto& pc = *ctx.phaseConstraints;

    if (!pc.hasActiveConstraints()) {
        return false;
    }

    // Count total number of constrained phases to validate against phase rule
    int nConstrainedSoln = 0;
    int nConstrainedCond = 0;
    for (const auto& c : pc.solnPhaseConstraints) {
        if (c.mode == PhaseConstraintMode::Fixed) {
            nConstrainedSoln++;
        }
    }
    for (const auto& c : pc.condPhaseConstraints) {
        if (c.mode == PhaseConstraintMode::Fixed) {
            nConstrainedCond++;
        }
    }

    // Count active elements (those with non-zero input mass)
    // The Gibbs phase rule limit is based on active elements, not total database elements
    int nActiveElements = 0;
    for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
        if (thermo.dMolesElement(j) > 0.0) {
            nActiveElements++;
        }
    }

    // Validate: total phases cannot exceed number of active elements (Gibbs phase rule)
    // Using nElements for array bounds, but nActiveElements for thermodynamic validity
    int nTotalConstrained = nConstrainedSoln + nConstrainedCond;
    if (nTotalConstrained > nActiveElements) {
        // Too many constrained phases for this system's active elements
        return false;
    }
    // Also check array bounds (should be redundant but safety check)
    if (nTotalConstrained > thermo.nElements) {
        return false;
    }

    // Clear existing assemblage
    thermo.nSolnPhases = 0;
    thermo.nConPhases = 0;
    for (int i = 0; i < thermo.nElements; ++i) {
        thermo.iAssemblage(i) = 0;
        thermo.dMolesPhase(i) = 0.0;
    }
    std::fill(gem.lSolnPhases.begin(), gem.lSolnPhases.end(), false);

    // Clear species moles and mole fractions to avoid stale values from
    // previously active phases contributing to mass balance calculations
    for (int i = 0; i < thermo.nSpecies; ++i) {
        thermo.dMolesSpecies(i) = 0.0;
        thermo.dMolFraction(i) = 0.0;
    }

    // Compute total element moles for initial estimates
    double totalElementMoles = 0.0;
    for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
        totalElementMoles += thermo.dMolesElement(j);
    }
    if (totalElementMoles <= 0.0) {
        return false;
    }

    // Add constrained solution phases to assemblage
    for (int iPhase = 0; iPhase < static_cast<int>(pc.solnPhaseConstraints.size()); ++iPhase) {
        auto& c = pc.solnPhaseConstraints[iPhase];
        if (c.mode != PhaseConstraintMode::Fixed) {
            continue;
        }

        // Add this solution phase
        thermo.nSolnPhases++;
        int newIdx = thermo.nElements - thermo.nSolnPhases;

        // Safety check: ensure index is valid (should be guaranteed by earlier validation)
        if (newIdx < thermo.nConPhases || newIdx < 0) {
            // Would overlap with condensed phases or go negative
            thermo.nSolnPhases--;
            return false;
        }

        thermo.iAssemblage(newIdx) = -(iPhase + 1);
        gem.lSolnPhases[iPhase] = true;

        // Estimate initial phase moles based on target fraction
        // phase_element_moles = target_fraction * total_element_moles
        // This is approximate; the solver will adjust
        double targetElementMoles = c.targetFraction * totalElementMoles;

        // Get species range for this phase
        int iFirst = (iPhase > 0) ? thermo.nSpeciesPhase(iPhase) : 0;
        int iLast = thermo.nSpeciesPhase(iPhase + 1);

        // Initialize mole fractions (equal distribution among feasible species)
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

        // Initialize mole fractions
        for (int i = iFirst; i < iLast; ++i) {
            thermo.dMolFraction(i) = 1e-100;
        }
        if (!feasibleSpecies.empty()) {
            double x = 1.0 / feasibleSpecies.size();
            for (int i : feasibleSpecies) {
                thermo.dMolFraction(i) = x;
            }
        }

        // Estimate phase moles from element content and stoichiometry
        // Use average stoichiometry over feasible species
        double avgStoich = 0.0;
        for (int i : feasibleSpecies) {
            for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
                avgStoich += thermo.dMolFraction(i) * thermo.dStoichSpecies(i, j);
            }
        }
        double phaseMoles = (avgStoich > 1e-10) ? targetElementMoles / avgStoich : 0.1;
        // Allow zero phase moles for zero target fractions (phase field use case)
        thermo.dMolesPhase(newIdx) = std::max(0.0, phaseMoles);

        // Update species moles
        for (int i = iFirst; i < iLast; ++i) {
            thermo.dMolesSpecies(i) = thermo.dMolesPhase(newIdx) * thermo.dMolFraction(i);
        }
    }

    // Add constrained pure condensed phases to assemblage
    int nSolnSpecies = thermo.nSpeciesPhase(thermo.nSolnPhasesSys);
    for (int iCond = 0; iCond < static_cast<int>(pc.condPhaseConstraints.size()); ++iCond) {
        auto& c = pc.condPhaseConstraints[iCond];
        if (c.mode != PhaseConstraintMode::Fixed) {
            continue;
        }

        int speciesIdx = nSolnSpecies + iCond;
        if (speciesIdx >= thermo.nSpecies) {
            continue;
        }

        // Safety check: ensure condensed phase index doesn't overlap with solution phases
        // Solution phases occupy indices [nElements - nSolnPhases, nElements)
        // Condensed phases occupy indices [0, nConPhases)
        if (thermo.nConPhases >= thermo.nElements - thermo.nSolnPhases) {
            return false;
        }

        // Add this pure condensed phase
        thermo.iAssemblage(thermo.nConPhases) = speciesIdx + 1;  // 1-based index

        // Estimate moles from target fraction
        double targetElementMoles = c.targetFraction * totalElementMoles;
        double stoichSum = 0.0;
        for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
            stoichSum += thermo.dStoichSpecies(speciesIdx, j);
        }
        double phaseMoles = (stoichSum > 1e-10) ? targetElementMoles / stoichSum : 0.1;
        // Allow zero phase moles for zero target fractions (phase field use case)
        thermo.dMolesPhase(thermo.nConPhases) = std::max(0.0, phaseMoles);
        thermo.dMolesSpecies(speciesIdx) = thermo.dMolesPhase(thermo.nConPhases);

        thermo.nConPhases++;
    }

    // Initialize element potentials if not already set
    // Use a simple estimate based on pure species Gibbs energies
    double R = Constants::kIdealGasConstant;
    double T = ctx.io->dTemperature;

    for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
        if (thermo.dMolesElement(j) <= 0) {
            thermo.dElementPotential(j) = 0.0;
            continue;
        }

        // Find a species containing this element to estimate potential
        double potentialSum = 0.0;
        int count = 0;
        for (int iPhase = 0; iPhase < thermo.nSolnPhases; ++iPhase) {
            int assembIdx = thermo.nElements - thermo.nSolnPhases + iPhase;
            int phaseIdx = -thermo.iAssemblage(assembIdx) - 1;
            if (phaseIdx < 0) continue;

            int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
            int iLast = thermo.nSpeciesPhase(phaseIdx + 1);

            for (int i = iFirst; i < iLast; ++i) {
                if (thermo.dStoichSpecies(i, j) > 0 && thermo.dMolFraction(i) > 1e-10) {
                    // λ_j ≈ (μ_std + ln(x)) * p / a_ij
                    double muStd = thermo.dStdGibbsEnergy(i) / (R * T);
                    double mu = muStd + std::log(thermo.dMolFraction(i));
                    double lambda = mu * thermo.iParticlesPerMole(i) / thermo.dStoichSpecies(i, j);
                    potentialSum += lambda;
                    count++;
                }
            }
        }

        if (count > 0) {
            thermo.dElementPotential(j) = potentialSum / count;
        }
    }

    return (thermo.nSolnPhases + thermo.nConPhases) > 0;
}

void ConstrainedGEM::updateConstrainedPhaseMoles(ThermoContext& ctx) {
    auto& thermo = *ctx.thermo;
    auto& io = *ctx.io;
    auto& pc = *ctx.phaseConstraints;

    double R = Constants::kIdealGasConstant;
    double T = io.dTemperature;

    // Compute total element moles
    double totalElementMoles = 0.0;
    for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
        totalElementMoles += thermo.dMolesElement(j);
    }
    if (totalElementMoles <= 0.0) return;

    // Update phase moles to match target fractions
    // Use damped iteration to avoid oscillation
    double damping = 0.3;

    for (int iPhase = 0; iPhase < thermo.nSolnPhases; ++iPhase) {
        int assembIdx = thermo.nElements - thermo.nSolnPhases + iPhase;
        int phaseIdx = -thermo.iAssemblage(assembIdx) - 1;
        if (phaseIdx < 0 || phaseIdx >= thermo.nSolnPhasesSys) continue;

        // Check if this phase is constrained
        if (phaseIdx >= static_cast<int>(pc.solnPhaseConstraints.size())) continue;
        auto& c = pc.solnPhaseConstraints[phaseIdx];
        if (c.mode != PhaseConstraintMode::Fixed) continue;

        // Get species range for this phase
        int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
        int iLast = thermo.nSpeciesPhase(phaseIdx + 1);

        // Compute current average stoichiometry (elements per mole of phase)
        double avgStoich = 0.0;
        for (int i = iFirst; i < iLast; ++i) {
            if (thermo.dMolFraction(i) > 0) {
                for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
                    avgStoich += thermo.dMolFraction(i) * thermo.dStoichSpecies(i, j);
                }
            }
        }

        // Target element moles for this phase
        double targetElementMoles = c.targetFraction * totalElementMoles;

        // Target phase moles = target element moles / average stoichiometry
        double minPhaseMoles = thermo.tolerances[kTolPhaseMoles];
        double targetPhaseMoles = (avgStoich > minPhaseMoles) ? targetElementMoles / avgStoich : 0.1;

        // Damped update toward target
        // Allow zero phase moles for zero target fractions (phase field use case)
        double currentMoles = thermo.dMolesPhase(assembIdx);
        double newMoles = currentMoles + damping * (targetPhaseMoles - currentMoles);
        thermo.dMolesPhase(assembIdx) = std::max(0.0, newMoles);

        // Update species moles
        double phaseMoles = thermo.dMolesPhase(assembIdx);
        for (int i = iFirst; i < iLast; ++i) {
            if (thermo.dMolFraction(i) > 0) {
                thermo.dMolesSpecies(i) = phaseMoles * thermo.dMolFraction(i);
            }
        }
    }

    // Update constrained pure condensed phase moles toward target fractions
    // Pure condensed phases have fixed stoichiometry, so we can only adjust total moles
    int nSolnSpecies = thermo.nSpeciesPhase(thermo.nSolnPhasesSys);
    for (int iCond = 0; iCond < thermo.nConPhases; ++iCond) {
        int speciesIdx = thermo.iAssemblage(iCond) - 1;  // 1-based to 0-based
        if (speciesIdx < 0 || speciesIdx >= thermo.nSpecies) continue;

        // Map from species index to condensed phase constraint index
        int condIdx = speciesIdx - nSolnSpecies;
        if (condIdx < 0 || condIdx >= static_cast<int>(pc.condPhaseConstraints.size())) continue;

        auto& c = pc.condPhaseConstraints[condIdx];
        if (c.mode != PhaseConstraintMode::Fixed) continue;

        // Compute stoichiometry sum for this species (elements per mole)
        double stoichSum = 0.0;
        for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
            stoichSum += thermo.dStoichSpecies(speciesIdx, j);
        }

        // Target element moles for this phase
        double targetElementMoles = c.targetFraction * totalElementMoles;

        // Target phase moles = target element moles / stoichiometry sum
        double minPhaseMoles = thermo.tolerances[kTolPhaseMoles];
        double targetPhaseMoles = (stoichSum > minPhaseMoles) ? targetElementMoles / stoichSum : 0.1;

        // Damped update toward target
        // Allow zero phase moles for zero target fractions (phase field use case)
        double currentMoles = thermo.dMolesPhase(iCond);
        double newMoles = currentMoles + damping * (targetPhaseMoles - currentMoles);
        thermo.dMolesPhase(iCond) = std::max(0.0, newMoles);

        // Update species moles (for pure condensed, species moles = phase moles)
        thermo.dMolesSpecies(speciesIdx) = thermo.dMolesPhase(iCond);
    }

    // Update mole fractions based on element potentials (for ideal mixing phases)
    // At equilibrium: x_i ∝ exp(μ* - μ_std) where μ* = Σ_j λ_j * a_{ij} / p_i
    for (int iPhase = 0; iPhase < thermo.nSolnPhases; ++iPhase) {
        int assembIdx = thermo.nElements - thermo.nSolnPhases + iPhase;
        int phaseIdx = -thermo.iAssemblage(assembIdx) - 1;
        if (phaseIdx < 0 || phaseIdx >= thermo.nSolnPhasesSys) continue;

        // Only update mole fractions for ideal mixing phases
        if (phaseIdx >= static_cast<int>(thermo.iSolnPhaseType.size()) ||
            thermo.iSolnPhaseType[phaseIdx] != Constants::PhaseType::IDMX) {
            continue;
        }

        int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
        int iLast = thermo.nSpeciesPhase(phaseIdx + 1);

        // Compute new mole fractions from element potentials
        double sumExp = 0.0;
        std::vector<double> expVal(iLast - iFirst, 0.0);

        for (int i = iFirst; i < iLast; ++i) {
            // Check feasibility
            bool feasible = true;
            for (int j = 0; j < thermo.nElements - thermo.nChargedConstraints; ++j) {
                if (thermo.dStoichSpecies(i, j) > 0.0 && thermo.dMolesElement(j) <= 0.0) {
                    feasible = false;
                    break;
                }
            }
            if (!feasible) continue;

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

        // Update mole fractions with damping
        if (sumExp > 1e-300) {
            double moleFracDamping = 0.5;
            for (int i = iFirst; i < iLast; ++i) {
                double xNew = expVal[i - iFirst] / sumExp;
                double xOld = thermo.dMolFraction(i);
                thermo.dMolFraction(i) = xOld + moleFracDamping * (xNew - xOld);
            }

            // Normalize to ensure sum = 1
            double sum = 0.0;
            for (int i = iFirst; i < iLast; ++i) {
                sum += thermo.dMolFraction(i);
            }
            if (sum > 1e-300) {
                for (int i = iFirst; i < iLast; ++i) {
                    thermo.dMolFraction(i) /= sum;
                }
            }
        }

        // Update species moles after mole fraction update
        double phaseMoles = thermo.dMolesPhase(assembIdx);
        for (int i = iFirst; i < iLast; ++i) {
            if (thermo.dMolFraction(i) > 0) {
                thermo.dMolesSpecies(i) = phaseMoles * thermo.dMolFraction(i);
            }
        }
    }

    // Update element potentials from chemical equilibrium
    // For species with single element: λ_j = μ_i * p_i / a_{ij}
    for (int iPhase = 0; iPhase < thermo.nSolnPhases; ++iPhase) {
        int assembIdx = thermo.nElements - thermo.nSolnPhases + iPhase;
        int phaseIdx = -thermo.iAssemblage(assembIdx) - 1;
        if (phaseIdx < 0 || phaseIdx >= thermo.nSolnPhasesSys) continue;

        int iFirst = (phaseIdx > 0) ? thermo.nSpeciesPhase(phaseIdx) : 0;
        int iLast = thermo.nSpeciesPhase(phaseIdx + 1);

        for (int i = iFirst; i < iLast; ++i) {
            if (thermo.dMolFraction(i) < 1e-20) continue;

            // Check if this is a pure end member (only one element)
            int singleElement = -1;
            int nElements = 0;
            for (int j = 0; j < thermo.nElements; ++j) {
                if (thermo.dStoichSpecies(i, j) > 0) {
                    singleElement = j;
                    nElements++;
                }
            }

            // For pure end member species, λ_j = μ_i * p_i / a_{ij}
            if (nElements == 1 && singleElement >= 0) {
                double muStd = thermo.dStdGibbsEnergy(i) / (R * T);
                double mu = muStd + std::log(thermo.dMolFraction(i));
                double stoich = thermo.dStoichSpecies(i, singleElement);
                double particles = thermo.iParticlesPerMole(i);
                thermo.dElementPotential(singleElement) = mu * particles / stoich;
            }
        }
    }
}

} // namespace Thermochimica
