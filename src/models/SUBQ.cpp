/// @file SUBQ.cpp
/// @brief SUBQ model variant of Modified Quasichemical Model
/// @details SUBQ is similar to SUBG but with different coordination number treatment

#include "thermochimica/ThermoContext.hpp"

namespace Thermochimica {

// Forward declaration of SUBG function
void computeExcessGibbsSUBG(ThermoContext& ctx, int phaseIndex);

/// @brief Compute partial molar excess Gibbs energy for SUBQ phases
/// @details SUBQ is a variant of the Modified Quasichemical Model (MQM) that uses
/// the same core algorithm as SUBG. The difference is in how coordination numbers
/// are specified in the database file.
///
/// @param ctx The thermochimica context
/// @param phaseIndex Absolute index of the solution phase
void computeExcessGibbsSUBQ(ThermoContext& ctx, int phaseIndex) {
    // SUBQ uses the same computational approach as SUBG
    // The difference is in the database format, not the calculation
    computeExcessGibbsSUBG(ctx, phaseIndex);
}

/// @brief Wrapper for backwards compatibility
void compExcessGibbsEnergySUBQ(ThermoContext& ctx, int iSolnIndex) {
    computeExcessGibbsSUBQ(ctx, iSolnIndex);
}

} // namespace Thermochimica
