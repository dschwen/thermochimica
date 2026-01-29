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

} // namespace Thermochimica
