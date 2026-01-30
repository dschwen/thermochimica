#include "thermochimica/solver/GEMSolver.hpp"
#include "thermochimica/context/PhaseConstraints.hpp"
#include "thermochimica/util/Constants.hpp"
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

    // Validate: total phases cannot exceed number of elements (Gibbs phase rule)
    // If exceeded, writes to iAssemblage/dMolesPhase would go out of bounds
    if (nConstrainedSoln + nConstrainedCond > thermo.nElements) {
        // Too many constrained phases for this system
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
        thermo.dMolesPhase(newIdx) = std::max(1e-10, phaseMoles);

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
        thermo.dMolesPhase(thermo.nConPhases) = std::max(1e-10, phaseMoles);
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

} // namespace Thermochimica
