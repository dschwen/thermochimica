#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>

namespace Thermochimica {

/// Common Tangent Zone state - converted from ModuleCTZ.f90
/// Data for CTZ interpolation calculations
struct CTZState {
    // Scalar integers
    int nMaxAssemblages = 0;    ///< Maximum assemblages stored
    int nMaxElements = 0;       ///< Maximum elements
    int nAssemblages = 0;       ///< Current assemblage count

    // Scalar reals
    double maxNorm = 0.0;       ///< Maximum norm tolerance
    double tRange = 0.0;        ///< Temperature range

    // Initialization flag
    bool lCtzInit = false;      ///< CTZ initialized flag

    // 2D integer array [assemblage x phases]
    Eigen::MatrixXi assemblageHistory;

    // 2D string array [assemblage x elements]
    std::vector<std::vector<std::string>> elementHistory;

    // 2D real array [assemblage x 2] for temperature limits
    Eigen::MatrixXd assemblageTlimits;

    // 4D real array stored as nested structure
    // [assemblage][phase][species][element] -> stoichiometry
    std::vector<std::vector<Eigen::MatrixXd>> stoichHistory;

    // Constructor
    CTZState() = default;

    /// Allocate CTZ storage
    void allocate(int maxAssemblages, int maxElements, int maxPhases, int maxSpecies);

    /// Initialize CTZ module
    void init(int maxAssemblages, int maxElements);

    /// Reset CTZ state
    void reset();

    /// Add an assemblage to history
    void addAssemblage(const Eigen::VectorXi& assemblage,
                       const std::vector<std::string>& elements,
                       double tMin, double tMax);
};

} // namespace Thermochimica
