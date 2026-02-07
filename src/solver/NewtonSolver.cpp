/// @file NewtonSolver.cpp
/// @brief Implementation of Newton direction solver

#include "thermochimica/solver/NewtonSolver.hpp"
#include "thermochimica/context/ThermoState.hpp"
#include "thermochimica/context/GEMState.hpp"
#include "thermochimica/context/ThermoIO.hpp"
#include "thermochimica/context/PhaseConstraints.hpp"
#include "thermochimica/solver/GEMSolver.hpp"
#include "thermochimica/util/Constants.hpp"
#include <Eigen/LU>
#include <cmath>
#include <iostream>

namespace Thermochimica {

NewtonSolver::NewtonSolver(ThermoIO& io, PhaseConstraints* constraints)
    : io_(io), constraints_(constraints) {
}

int NewtonSolver::computeDirection(ThermoState& state,
                                   GEMState& gemState,
                                   Eigen::VectorXd& direction) {
    // Build list of active elements (with non-zero input mass)
    std::vector<int> activeElements;
    buildActiveElements(state, activeElements);
    int nActiveEl = static_cast<int>(activeElements.size());

    // Determine reduced system size (active elements + phases)
    int nVar = nActiveEl + state.nConPhases + state.nSolnPhases;

    // Construct Hessian and RHS for reduced system
    Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(nVar, nVar);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(nVar);

    constructReducedSystem(state, hessian, rhs, activeElements);

    // Add phase constraint penalty contributions to Newton system
    if (constraints_ != nullptr && constraints_->hasActiveConstraints()) {
        addConstraintContributions(state, gemState, hessian, rhs, activeElements);
    }

    // Solve using Eigen's PartialPivLU (equivalent to LAPACK DGESV)
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(hessian);

    // Check for singularity
    if (std::abs(lu.determinant()) < 1e-20) {
        handleSingular(hessian, rhs);
        lu.compute(hessian);
    }

    // Solve the reduced system
    Eigen::VectorXd updateReduced = lu.solve(rhs);

    // Map back to full update vector
    mapToFullDirection(updateReduced, direction, state, activeElements);

    // Check for NaN in solution
    int nVarFull = state.nElements + state.nConPhases + state.nSolnPhases;
    for (int i = 0; i < nVarFull; ++i) {
        if (std::isnan(direction(i)) || std::isinf(direction(i))) {
            gemState.lRevertSystem = true;
            return 1;
        }
    }

    return 0;
}

void NewtonSolver::handleSingular(Eigen::MatrixXd& hessian, Eigen::VectorXd& rhs) {
    // For singular matrices, add small diagonal perturbation
    int n = hessian.rows();
    for (int i = 0; i < n; ++i) {
        if (std::abs(hessian(i, i)) < 1e-20) {
            hessian(i, i) = 1.0;
            rhs(i) = 0.0;
        }
    }
}

void NewtonSolver::buildActiveElements(const ThermoState& state,
                                       std::vector<int>& activeElements) const {
    activeElements.clear();
    for (int j = 0; j < state.nElements - state.nChargedConstraints; ++j) {
        if (state.dMolesElement(j) > 0.0) {
            activeElements.push_back(j);
        }
    }
}

void NewtonSolver::constructReducedSystem(const ThermoState& state,
                                          Eigen::MatrixXd& hessian,
                                          Eigen::VectorXd& rhs,
                                          const std::vector<int>& activeElements) const {
    int nActiveEl = static_cast<int>(activeElements.size());
    int nConPhases = state.nConPhases;
    int nSolnPhases = state.nSolnPhases;
    int nElements = state.nElements;
    double R = Constants::kIdealGasConstant;
    double T = io_.dTemperature;

    // Build Hessian matrix blocks (using reduced element indexing)
    // Block 1: Element-element block (from solution phase species)
    for (int iPhase = 0; iPhase < nSolnPhases; ++iPhase) {
        int phaseIdx = -state.iAssemblage(nElements - nSolnPhases + iPhase) - 1;
        if (phaseIdx < 0) continue;

        int iFirst = (phaseIdx > 0) ? state.nSpeciesPhase(phaseIdx) : 0;
        int iLast = state.nSpeciesPhase(phaseIdx + 1);

        for (int iSpecies = iFirst; iSpecies < iLast; ++iSpecies) {
            if (state.dMolFraction(iSpecies) <= 0) continue;

            double molesSpec = state.dMolesSpecies(iSpecies);
            int particles = state.iParticlesPerMole(iSpecies);
            double factor = molesSpec / (particles * particles);

            for (int jRed = 0; jRed < nActiveEl; ++jRed) {
                int jFull = activeElements[jRed];
                for (int iRed = jRed; iRed < nActiveEl; ++iRed) {
                    int iFull = activeElements[iRed];
                    double contrib = state.dStoichSpecies(iSpecies, iFull) *
                                    state.dStoichSpecies(iSpecies, jFull) * factor;
                    hessian(iRed, jRed) += contrib;
                    if (iRed != jRed) {
                        hessian(jRed, iRed) += contrib;  // Symmetry
                    }
                }
            }
        }
    }

    // Block 2: Element-Solution phase block
    for (int iPhase = 0; iPhase < nSolnPhases; ++iPhase) {
        int col = nActiveEl + iPhase;
        int phaseIdx = -state.iAssemblage(nElements - nSolnPhases + iPhase) - 1;
        if (phaseIdx < 0) continue;

        int iFirst = (phaseIdx > 0) ? state.nSpeciesPhase(phaseIdx) : 0;
        int iLast = state.nSpeciesPhase(phaseIdx + 1);

        for (int jRed = 0; jRed < nActiveEl; ++jRed) {
            int jFull = activeElements[jRed];
            double sum = 0.0;
            for (int iSpecies = iFirst; iSpecies < iLast; ++iSpecies) {
                if (state.dMolFraction(iSpecies) > 0) {
                    sum += state.dStoichSpecies(iSpecies, jFull) *
                           state.dMolFraction(iSpecies) /
                           state.iParticlesPerMole(iSpecies);
                }
            }
            hessian(jRed, col) = sum;
            hessian(col, jRed) = sum;  // Symmetry
        }
    }

    // Block 3: Element-Condensed phase block
    for (int iPhase = 0; iPhase < nConPhases; ++iPhase) {
        int col = nActiveEl + nSolnPhases + iPhase;
        int speciesIdx = state.iAssemblage(iPhase) - 1;
        if (speciesIdx < 0 || speciesIdx >= state.nSpecies) continue;

        for (int jRed = 0; jRed < nActiveEl; ++jRed) {
            int jFull = activeElements[jRed];
            double val = state.dStoichSpecies(speciesIdx, jFull) /
                        state.iParticlesPerMole(speciesIdx);
            hessian(jRed, col) = val;
            hessian(col, jRed) = val;
        }
    }

    // Construct RHS as residuals (negative gradient for Newton)
    // Newton solves: H * Δx = -r, where r is the residual

    // Element mass balance residuals (for active elements only)
    // Residual = target - current = b_j - ∑_i (n_i * a_ij / p_i)
    for (int jRed = 0; jRed < nActiveEl; ++jRed) {
        int jFull = activeElements[jRed];
        double target = state.dMolesElement(jFull);
        double current = 0.0;

        // Sum contribution from solution phase species
        for (int iPhase = 0; iPhase < nSolnPhases; ++iPhase) {
            int phaseIdx = -state.iAssemblage(nElements - nSolnPhases + iPhase) - 1;
            if (phaseIdx < 0) continue;

            int iFirst = (phaseIdx > 0) ? state.nSpeciesPhase(phaseIdx) : 0;
            int iLast = state.nSpeciesPhase(phaseIdx + 1);

            for (int iSpecies = iFirst; iSpecies < iLast; ++iSpecies) {
                if (state.dMolFraction(iSpecies) > 0) {
                    current += state.dMolesSpecies(iSpecies) *
                               state.dStoichSpecies(iSpecies, jFull) /
                               state.iParticlesPerMole(iSpecies);
                }
            }
        }

        // Sum contribution from pure condensed phase species
        for (int iPhase = 0; iPhase < nConPhases; ++iPhase) {
            int speciesIdx = state.iAssemblage(iPhase) - 1;
            if (speciesIdx < 0 || speciesIdx >= state.nSpecies) continue;

            current += state.dMolesPhase(iPhase) *
                       state.dStoichSpecies(speciesIdx, jFull) /
                       state.iParticlesPerMole(speciesIdx);
        }

        // RHS = negative residual for Newton method
        rhs(jRed) = -(target - current);
    }

    // Solution phase Gibbs energy residuals
    // At equilibrium: dGibbsSolnPhase = 0 (normalized driving force)
    // Residual = 0 - current = -dGibbsSolnPhase
    for (int iPhase = 0; iPhase < nSolnPhases; ++iPhase) {
        int row = nActiveEl + iPhase;
        int phaseIdx = -state.iAssemblage(nElements - nSolnPhases + iPhase) - 1;
        if (phaseIdx < 0) continue;

        // RHS = negative residual = -(-dGibbsSolnPhase) = dGibbsSolnPhase
        rhs(row) = state.dGibbsSolnPhase(phaseIdx);
    }

    // Pure condensed phase equilibrium residuals
    // At equilibrium: G°/RT = ∑_j λ_j * a_ij / p_i
    // Residual = G°/RT - ∑_j λ_j * a_ij / p_i
    for (int iPhase = 0; iPhase < nConPhases; ++iPhase) {
        int row = nActiveEl + nSolnPhases + iPhase;
        int speciesIdx = state.iAssemblage(iPhase) - 1;
        if (speciesIdx < 0 || speciesIdx >= state.nSpecies) continue;

        double gRT = state.dStdGibbsEnergy(speciesIdx) / (R * T);

        // Compute current ∑_j λ_j * a_ij / p_i
        double lambdaSum = 0.0;
        for (int jRed = 0; jRed < nActiveEl; ++jRed) {
            int jFull = activeElements[jRed];
            lambdaSum += state.dElementPotential(jFull) *
                         state.dStoichSpecies(speciesIdx, jFull);
        }
        lambdaSum /= state.iParticlesPerMole(speciesIdx);

        // RHS = negative residual for Newton method
        rhs(row) = -(gRT - lambdaSum);
    }
}

void NewtonSolver::addConstraintContributions(const ThermoState& state,
                                              const GEMState& gemState,
                                              Eigen::MatrixXd& hessian,
                                              Eigen::VectorXd& rhs,
                                              const std::vector<int>& activeElements) const {
    // Create temporary context for bridge to legacy code
    ThermoContext ctx;
    ThermoState* statePtr = const_cast<ThermoState*>(&state);
    GEMState* gemPtr = const_cast<GEMState*>(&gemState);
    ctx.thermo.reset(statePtr);
    ctx.gem.reset(gemPtr);
    ctx.phaseConstraints.reset(const_cast<PhaseConstraints*>(constraints_));

    // Call legacy functions to add constraint contributions
    ConstrainedGEM::addConstraintGradients(ctx, rhs, activeElements);
    ConstrainedGEM::addConstraintHessian(ctx, hessian, activeElements);

    // Release pointers so they won't be deleted
    ctx.thermo.release();
    ctx.gem.release();
    ctx.phaseConstraints.release();
}

void NewtonSolver::mapToFullDirection(const Eigen::VectorXd& updateReduced,
                                      Eigen::VectorXd& direction,
                                      const ThermoState& state,
                                      const std::vector<int>& activeElements) const {
    int nActiveEl = static_cast<int>(activeElements.size());

    // Initialize full direction to zero
    int nVarFull = state.nElements + state.nConPhases + state.nSolnPhases;
    direction.resize(nVarFull);
    direction.setZero();

    // Map active element updates
    for (int i = 0; i < nActiveEl; ++i) {
        direction(activeElements[i]) = updateReduced(i);
    }

    // Solution phase updates (map from reduced to full indexing)
    for (int i = 0; i < state.nSolnPhases; ++i) {
        direction(state.nElements + i) = updateReduced(nActiveEl + i);
    }

    // Condensed phase updates
    for (int i = 0; i < state.nConPhases; ++i) {
        direction(state.nElements + state.nSolnPhases + i) =
            updateReduced(nActiveEl + state.nSolnPhases + i);
    }
}

} // namespace Thermochimica
