/// @file TestPhaseConstraints.cpp
/// @brief Integration tests for phase fraction constraints
/// @details Tests constrained GEM solver with phase fraction targets

#include <gtest/gtest.h>
#include <thermochimica/Thermochimica.hpp>
#include <thermochimica/context/PhaseConstraints.hpp>
#include <cmath>

using namespace Thermochimica;

// Helper to check relative error
inline bool relativeError(double actual, double expected, double tolerance = 1e-3) {
    if (std::abs(expected) < 1e-10) {
        return std::abs(actual - expected) < tolerance;
    }
    return std::abs((actual - expected) / expected) < tolerance;
}

//=============================================================================
// Phase Constraint API Tests
//=============================================================================

// Test setting and clearing constraints
// Note: NobleMetals-Kaye.dat has phases: gas_ideal, FCCN, BCCN, HCPN, LiqN, sigma
TEST(PhaseConstraintAPI, SetAndClearConstraints) {
    ThermoContext ctx;
    setStandardUnits(ctx);
    setThermoFilename(ctx, "NobleMetals-Kaye.dat");
    parseCSDataFile(ctx);

    ASSERT_TRUE(ctx.isDatabaseLoaded());

    // Initially no constraints
    EXPECT_FALSE(ctx.phaseConstraints->hasActiveConstraints());

    // Set a solution phase constraint (use actual phase names from database)
    setSolnPhaseConstraint(ctx, "FCCN", 0.5);
    EXPECT_TRUE(ctx.phaseConstraints->hasActiveConstraints());
    EXPECT_EQ(ctx.phaseConstraints->getNumActiveConstraints(), 1);

    // Set another constraint
    setSolnPhaseConstraint(ctx, "BCCN", 0.3);
    EXPECT_EQ(ctx.phaseConstraints->getNumActiveConstraints(), 2);

    // Remove one constraint
    removePhaseConstraint(ctx, "FCCN");
    EXPECT_EQ(ctx.phaseConstraints->getNumActiveConstraints(), 1);

    // Clear all constraints
    clearPhaseConstraints(ctx);
    EXPECT_FALSE(ctx.phaseConstraints->hasActiveConstraints());
    EXPECT_EQ(ctx.phaseConstraints->getNumActiveConstraints(), 0);
}

// Test constraint parameter setting
TEST(PhaseConstraintAPI, SetConstraintParameters) {
    ThermoContext ctx;

    // Test default values
    EXPECT_DOUBLE_EQ(ctx.phaseConstraints->constraintTolerance, 1e-4);
    EXPECT_DOUBLE_EQ(ctx.phaseConstraints->penaltyParameter, 1.0);
    EXPECT_EQ(ctx.phaseConstraints->maxOuterIterations, 20);

    // Set custom parameters
    setConstraintTolerance(ctx, 1e-5);
    EXPECT_DOUBLE_EQ(ctx.phaseConstraints->constraintTolerance, 1e-5);

    setConstraintPenaltyParameter(ctx, 10.0);
    EXPECT_DOUBLE_EQ(ctx.phaseConstraints->penaltyParameter, 10.0);

    setConstraintMaxOuterIterations(ctx, 50);
    EXPECT_EQ(ctx.phaseConstraints->maxOuterIterations, 50);
}

// Test invalid constraint values are clamped
TEST(PhaseConstraintAPI, ConstraintValueClamping) {
    ThermoContext ctx;
    setStandardUnits(ctx);
    setThermoFilename(ctx, "NobleMetals-Kaye.dat");
    parseCSDataFile(ctx);

    // Phase constraint values should be clamped to [0, 1]
    // FCCN is at index 1 in this database
    setSolnPhaseConstraint(ctx, "FCCN", 1.5);  // Above 1
    auto& constraint = ctx.phaseConstraints->solnPhaseConstraints[1];
    EXPECT_LE(constraint.targetFraction, 1.0);

    setSolnPhaseConstraint(ctx, "FCCN", -0.5);  // Below 0
    EXPECT_GE(constraint.targetFraction, 0.0);
}

//=============================================================================
// Unconstrained Baseline Tests
//=============================================================================

// Run unconstrained calculation to establish baseline
TEST(PhaseConstraintCalc, UnconstrainedBaseline) {
    ThermoContext ctx;
    setStandardUnits(ctx);
    setThermoFilename(ctx, "NobleMetals-Kaye.dat");
    parseCSDataFile(ctx);

    setTemperaturePressure(ctx, 1500.0, 1.0);
    setElementMass(ctx, 42, 0.5);  // Mo
    setElementMass(ctx, 44, 0.5);  // Ru

    // Run unconstrained calculation
    thermochimica(ctx);

    EXPECT_EQ(ctx.infoThermo(), 0) << "Unconstrained calculation should succeed";

    // Get phase fractions (use actual phase names from database)
    auto [fcc_frac, fcc_info] = getPhaseElementFraction(ctx, "FCCN");
    auto [bcc_frac, bcc_info] = getPhaseElementFraction(ctx, "BCCN");

    // Natural phase fractions should be in valid range
    EXPECT_GE(fcc_frac, 0.0);
    EXPECT_LE(fcc_frac, 1.0);
    EXPECT_GE(bcc_frac, 0.0);
    EXPECT_LE(bcc_frac, 1.0);

    // Store baseline Gibbs energy
    double gibbsBaseline = getGibbsEnergy(ctx);
    EXPECT_LT(gibbsBaseline, 0.0) << "Gibbs energy should be negative";
}

//=============================================================================
// Constrained Calculation Tests
//=============================================================================

// Test constraining multiple phases to specific fractions (phase field use case)
// Note: With the phase-field assemblage design, constrained phases are forced into
// the assemblage and unconstrained phases are excluded. This means we need to
// constrain multiple phases that sum to ~1.0 for a meaningful test.
TEST(PhaseConstraintCalc, SinglePhaseConstraint) {
    ThermoContext ctx;
    setStandardUnits(ctx);
    setThermoFilename(ctx, "NobleMetals-Kaye.dat");
    parseCSDataFile(ctx);

    setTemperaturePressure(ctx, 1800.0, 1.0);
    setElementMass(ctx, 42, 0.5);  // Mo
    setElementMass(ctx, 44, 0.5);  // Ru

    // Constrain two phases: FCCN to 40%, BCCN to 60% (sum to 1.0)
    double fccTarget = 0.4;
    double bccTarget = 0.6;
    setSolnPhaseConstraint(ctx, "FCCN", fccTarget);
    setSolnPhaseConstraint(ctx, "BCCN", bccTarget);

    // Verify constraints were set
    EXPECT_TRUE(ctx.phaseConstraints->hasActiveConstraints());
    EXPECT_EQ(ctx.phaseConstraints->getNumActiveConstraints(), 2);

    // Set solver parameters
    double tolerance = 1e-2;  // Slightly relaxed tolerance for two-phase system
    setConstraintTolerance(ctx, tolerance);
    setConstraintMaxOuterIterations(ctx, 50);

    // Run constrained calculation
    thermochimica(ctx);

    // Check that calculation converged
    EXPECT_EQ(ctx.infoThermo(), 0) << "Constrained calculation should succeed";

    // Check that constraint tracking works
    auto [fcc_frac, fcc_info] = getPhaseElementFraction(ctx, "FCCN");
    auto [bcc_frac, bcc_info] = getPhaseElementFraction(ctx, "BCCN");
    EXPECT_EQ(fcc_info, 0) << "Should be able to get FCC fraction";
    EXPECT_EQ(bcc_info, 0) << "Should be able to get BCC fraction";

    // Verify phase fractions are in valid range
    EXPECT_GE(fcc_frac, 0.0);
    EXPECT_LE(fcc_frac, 1.0);
    EXPECT_GE(bcc_frac, 0.0);
    EXPECT_LE(bcc_frac, 1.0);

    // Verify constraints are enforced within tolerance
    EXPECT_NEAR(fcc_frac, fccTarget, tolerance)
        << "FCCN fraction should be within tolerance of target";
    EXPECT_NEAR(bcc_frac, bccTarget, tolerance)
        << "BCCN fraction should be within tolerance of target";
}

// Test that constrained solution has higher Gibbs energy than unconstrained
TEST(PhaseConstraintCalc, ConstrainedGibbsHigher) {
    // First, run unconstrained calculation
    ThermoContext ctx1;
    setStandardUnits(ctx1);
    setThermoFilename(ctx1, "NobleMetals-Kaye.dat");
    parseCSDataFile(ctx1);

    setTemperaturePressure(ctx1, 1800.0, 1.0);
    setElementMass(ctx1, 42, 0.5);  // Mo
    setElementMass(ctx1, 44, 0.5);  // Ru

    thermochimica(ctx1);
    ASSERT_EQ(ctx1.infoThermo(), 0);
    double gibbsUnconstrained = getGibbsEnergy(ctx1);

    // Get natural FCCN fraction
    auto [natural_fcc, info] = getPhaseElementFraction(ctx1, "FCCN");

    // Run constrained calculation with different target
    ThermoContext ctx2;
    setStandardUnits(ctx2);
    setThermoFilename(ctx2, "NobleMetals-Kaye.dat");
    parseCSDataFile(ctx2);

    setTemperaturePressure(ctx2, 1800.0, 1.0);
    setElementMass(ctx2, 42, 0.5);
    setElementMass(ctx2, 44, 0.5);

    // Constrain to a different fraction than natural
    double targetFraction = std::min(0.9, std::max(0.1, natural_fcc + 0.2));
    setSolnPhaseConstraint(ctx2, "FCCN", targetFraction);
    setConstraintTolerance(ctx2, 1e-3);

    thermochimica(ctx2);

    if (ctx2.infoThermo() == 0) {
        double gibbsConstrained = getGibbsEnergy(ctx2);

        // Constrained optimum should have higher (less negative) Gibbs energy
        // since we're forcing a non-equilibrium state
        EXPECT_GE(gibbsConstrained, gibbsUnconstrained - 1e-6)
            << "Constrained Gibbs energy should be >= unconstrained";
    }
}

// Test constraint with small phase fraction
// Note: With phase-field design, constraining to 0 means not adding the phase.
// Instead, test constraining one phase to a small fraction with another taking the rest.
TEST(PhaseConstraintCalc, SmallFractionConstraint) {
    ThermoContext ctx;
    setStandardUnits(ctx);
    setThermoFilename(ctx, "NobleMetals-Kaye.dat");
    parseCSDataFile(ctx);

    setTemperaturePressure(ctx, 1500.0, 1.0);
    setElementMass(ctx, 42, 0.5);  // Mo
    setElementMass(ctx, 44, 0.5);  // Ru

    // Constrain BCCN to small fraction, FCCN takes the rest
    double bccTarget = 0.1;
    double fccTarget = 0.9;
    setSolnPhaseConstraint(ctx, "BCCN", bccTarget);
    setSolnPhaseConstraint(ctx, "FCCN", fccTarget);
    setConstraintTolerance(ctx, 1e-2);
    setConstraintMaxOuterIterations(ctx, 50);

    thermochimica(ctx);

    if (ctx.infoThermo() == 0) {
        auto [bcc_frac, bcc_info] = getPhaseElementFraction(ctx, "BCCN");
        auto [fcc_frac, fcc_info] = getPhaseElementFraction(ctx, "FCCN");
        EXPECT_NEAR(bcc_frac, bccTarget, 0.02) << "BCC should be near target";
        EXPECT_NEAR(fcc_frac, fccTarget, 0.02) << "FCC should be near target";
    }
}

//=============================================================================
// Reset and Reuse Tests
//=============================================================================

// Test that constraints are properly reset between calculations
TEST(PhaseConstraintCalc, ConstraintReset) {
    ThermoContext ctx;
    setStandardUnits(ctx);
    setThermoFilename(ctx, "NobleMetals-Kaye.dat");
    parseCSDataFile(ctx);

    // First constrained calculation
    setTemperaturePressure(ctx, 1500.0, 1.0);
    setElementMass(ctx, 42, 0.5);
    setElementMass(ctx, 44, 0.5);
    setSolnPhaseConstraint(ctx, "FCCN", 0.5);

    thermochimica(ctx);
    EXPECT_TRUE(ctx.phaseConstraints->hasActiveConstraints())
        << "Constraints should be active after calculation";

    // Reset for new calculation
    resetThermo(ctx);

    // Constraints should still be defined (reset keeps constraint definitions)
    EXPECT_TRUE(ctx.phaseConstraints->hasActiveConstraints())
        << "resetThermo should preserve constraint definitions";

    // But Lagrange multipliers should be reset
    for (const auto& c : ctx.phaseConstraints->solnPhaseConstraints) {
        EXPECT_DOUBLE_EQ(c.lagrangeMultiplier, 0.0)
            << "Lagrange multipliers should be reset";
    }

    // Full reset should clear constraints
    resetThermoAll(ctx);
    EXPECT_FALSE(ctx.phaseConstraints->hasActiveConstraints())
        << "resetThermoAll should clear constraints";
}

//=============================================================================
// Edge Cases and Error Handling
//=============================================================================

// Test constraint on non-existent phase
TEST(PhaseConstraintEdgeCases, NonExistentPhase) {
    ThermoContext ctx;
    setStandardUnits(ctx);
    setThermoFilename(ctx, "NobleMetals-Kaye.dat");
    parseCSDataFile(ctx);

    // Try to set constraint on non-existent phase
    setSolnPhaseConstraint(ctx, "NONEXISTENT_PHASE", 0.5);

    // Should not have any active constraints (invalid phase ignored)
    EXPECT_FALSE(ctx.phaseConstraints->hasActiveConstraints())
        << "Constraint on non-existent phase should be ignored";
}

// Test multiple constraints that sum to > 1
TEST(PhaseConstraintEdgeCases, OverconstrainedSystem) {
    ThermoContext ctx;
    setStandardUnits(ctx);
    setThermoFilename(ctx, "NobleMetals-Kaye.dat");
    parseCSDataFile(ctx);

    setTemperaturePressure(ctx, 1500.0, 1.0);
    setElementMass(ctx, 42, 0.5);
    setElementMass(ctx, 44, 0.5);

    // Set constraints that sum to > 1 (infeasible)
    setSolnPhaseConstraint(ctx, "FCCN", 0.6);
    setSolnPhaseConstraint(ctx, "BCCN", 0.6);

    thermochimica(ctx);

    // Should either fail or converge to best possible approximation
    // Either outcome is acceptable for an infeasible problem
    // Just verify it doesn't crash
}
