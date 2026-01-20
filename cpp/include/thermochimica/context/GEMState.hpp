#pragma once

#include <Eigen/Dense>
#include <vector>
#include "../util/Constants.hpp"

namespace Thermochimica {

/// GEM Solver state - converted from ModuleGEMSolver.f90
/// Stores state variables for Gibbs Energy Minimization solver iterations
struct GEMState {
    // Iteration tracking
    int iterLast = 0;               ///< Last phase assemblage change iteration
    int iterStep = 0;               ///< Current step
    int iterRevert = 0;             ///< Last reversion iteration
    int iterGlobal = 0;             ///< Current global iteration
    int iterLastCon = 0;            ///< Last condensed phase change
    int iterLastSoln = 0;           ///< Last solution phase change
    int iterSwap = 0;               ///< Swap counter
    int iterLastMiscGapCheck = 0;   ///< Last miscibility gap check

    // Phase tracking
    int iConPhaseLast = 0;          ///< Last condensed phase index
    int iSolnPhaseLast = 0;         ///< Last solution phase index
    int iSolnSwap = 0;              ///< Solution swap index
    int iPureConSwap = 0;           ///< Pure condensed swap index

    // Iteration history matrix [nElements x iterGlobalMax]
    Eigen::MatrixXi iterHistory;

    // Convergence metrics
    double dGEMFunctionNorm = 0.0;      ///< Current function norm
    double dGEMFunctionNormLast = 0.0;  ///< Previous function norm
    double dMaxSpeciesChange = 0.0;     ///< Maximum species change
    double dMinGibbs = 0.0;             ///< Minimum Gibbs energy

    // Working arrays
    Eigen::VectorXd dSumMolFractionSoln;   ///< Mole fraction sums per phase
    Eigen::VectorXd dMolesPhaseLast;       ///< Previous phase moles
    Eigen::VectorXd dUpdateVar;            ///< Direction/update vector
    Eigen::VectorXd dDrivingForceSoln;     ///< Driving forces per phase
    Eigen::VectorXd dPartialExcessGibbs;   ///< Partial excess Gibbs
    Eigen::VectorXd dPartialExcessGibbsLast; ///< Previous partial excess Gibbs

    // 2D working matrix
    Eigen::MatrixXd dEffStoichSolnPhase;   ///< Effective stoichiometry

    // Status flags
    bool lDebugMode = false;        ///< Debug output enabled
    bool lRevertSystem = false;     ///< System reversion needed
    bool lConverged = false;        ///< Convergence achieved

    // Phase status arrays
    std::vector<bool> lSolnPhases;  ///< Solution phase stability
    std::vector<bool> lMiscibility; ///< Miscibility gap flags

    // Constructor
    GEMState() = default;

    /// Allocate arrays based on system dimensions
    void allocate(int nElements, int nSolnPhasesSys, int nSpecies);

    /// Reset solver state for new calculation
    void reset();

    /// Initialize for a new iteration
    void initIteration();
};

} // namespace Thermochimica
