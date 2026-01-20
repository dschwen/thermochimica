#pragma once

#include <Eigen/Dense>
#include <vector>

namespace Thermochimica {

/// Subminimization state - converted from ModuleSubMin.f90
/// State for subminimization routine for solution phase composition
struct SubMinState {
    // Constants
    static constexpr double kSubMinTolerance = 1.0e-3;
    static constexpr double kMinMoleFraction = 1.0e-100;
    static constexpr double kTolEuclideanNorm = 1.0e-2;
    static constexpr double kTolDrivingForceChange = 1.0e-3;

    // Scalar integers
    int nVar = 0;                       ///< Number of variables in phase
    int iFirst = 0;                     ///< First species in phase
    int iLast = 0;                      ///< Last species in phase
    int iSolnPhaseIndexOther = 0;       ///< Reference phase index

    // Hessian structure
    Eigen::VectorXi iHessian;           ///< Hessian matrix structure indices

    // Scalar reals
    double dDrivingForce = 0.0;         ///< Current driving force
    double dDrivingForceLast = 0.0;     ///< Previous driving force
    double dSubMinFunctionNorm = 0.0;   ///< Function norm

    // Working arrays
    Eigen::VectorXd dChemicalPotentialStar;  ///< Chemical potentials from element potentials
    Eigen::VectorXd dRHS;                    ///< Working vector (functional/direction)
    Eigen::MatrixXd dHessian;                ///< Hessian matrix

    // Status
    bool lSubMinConverged = false;      ///< Convergence flag

    // Constructor
    SubMinState() = default;

    /// Allocate arrays for a phase with nSpecies species
    void allocate(int nSpecies);

    /// Reset for new subminimization
    void reset();

    /// Initialize for a specific phase
    void initForPhase(int first, int last);
};

} // namespace Thermochimica
