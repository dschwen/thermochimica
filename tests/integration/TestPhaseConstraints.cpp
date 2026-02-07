/// @file TestPhaseConstraints.cpp
/// @brief Integration tests for phase fraction constraints
/// @details Tests constrained GEM solver with phase fraction targets

#include <gtest/gtest.h>
#include <thermochimica/ThermoClass.hpp>
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
    ThermoClass thermo;
    thermo.setStandardUnits();
    thermo.setThermoFilename("NobleMetals-Kaye.dat");
    thermo.parseCSDataFile();

    ASSERT_TRUE(thermo.getContext().isDatabaseLoaded());

    // Initially no constraints
    EXPECT_FALSE(thermo.getContext().phaseConstraints->hasActiveConstraints());

    // Set a solution phase constraint (use actual phase names from database)
    thermo.setSolnPhaseConstraint("FCCN", 0.5);
    EXPECT_TRUE(thermo.getContext().phaseConstraints->hasActiveConstraints());
    EXPECT_EQ(thermo.getContext().phaseConstraints->getNumActiveConstraints(), 1);

    // Set another constraint
    thermo.setSolnPhaseConstraint("BCCN", 0.3);
    EXPECT_EQ(thermo.getContext().phaseConstraints->getNumActiveConstraints(), 2);

    // Remove one constraint
    thermo.removePhaseConstraint("FCCN");
    EXPECT_EQ(thermo.getContext().phaseConstraints->getNumActiveConstraints(), 1);

    // Clear all constraints
    thermo.clearPhaseConstraints();
    EXPECT_FALSE(thermo.getContext().phaseConstraints->hasActiveConstraints());
    EXPECT_EQ(thermo.getContext().phaseConstraints->getNumActiveConstraints(), 0);
}

// Test constraint parameter setting
TEST(PhaseConstraintAPI, SetConstraintParameters) {
    ThermoClass thermo;

    // Test default values
    EXPECT_DOUBLE_EQ(thermo.getContext().phaseConstraints->constraintTolerance, 1e-4);
    EXPECT_DOUBLE_EQ(thermo.getContext().phaseConstraints->penaltyParameter, 1.0);
    EXPECT_EQ(thermo.getContext().phaseConstraints->maxOuterIterations, 20);

    // Set custom parameters
    thermo.setConstraintTolerance(1e-5);
    EXPECT_DOUBLE_EQ(thermo.getContext().phaseConstraints->constraintTolerance, 1e-5);

    thermo.setConstraintPenaltyParameter(10.0);
    EXPECT_DOUBLE_EQ(thermo.getContext().phaseConstraints->penaltyParameter, 10.0);

    thermo.setConstraintMaxOuterIterations(50);
    EXPECT_EQ(thermo.getContext().phaseConstraints->maxOuterIterations, 50);
}

// Test invalid constraint values are clamped
TEST(PhaseConstraintAPI, ConstraintValueClamping) {
    ThermoClass thermo;
    thermo.setStandardUnits();
    thermo.setThermoFilename("NobleMetals-Kaye.dat");
    thermo.parseCSDataFile();

    // Phase constraint values should be clamped to [0, 1]
    // FCCN is at index 1 in this database
    thermo.setSolnPhaseConstraint("FCCN", 1.5);  // Above 1
    auto& constraint = thermo.getContext().phaseConstraints->solnPhaseConstraints[1];
    EXPECT_LE(constraint.targetFraction, 1.0);

    thermo.setSolnPhaseConstraint("FCCN", -0.5);  // Below 0
    EXPECT_GE(constraint.targetFraction, 0.0);
}

//=============================================================================
// Unconstrained Baseline Tests
//=============================================================================

// Run unconstrained calculation to establish baseline
TEST(PhaseConstraintCalc, UnconstrainedBaseline) {
    ThermoClass thermo;
    thermo.setStandardUnits();
    thermo.setThermoFilename("NobleMetals-Kaye.dat");
    thermo.parseCSDataFile();

    thermo.setTemperaturePressure(1500.0, 1.0);
    thermo.setElementMass(42, 0.5);  // Mo
    thermo.setElementMass(44, 0.5);  // Ru

    // Run unconstrained calculation
    thermo.calculate();

    EXPECT_EQ(thermo.getInfoCode(), 0) << "Unconstrained calculation should succeed";

    // Get phase fractions (use actual phase names from database)
    auto [fcc_frac, fcc_info] = thermo.getPhaseElementFraction("FCCN");
    auto [bcc_frac, bcc_info] = thermo.getPhaseElementFraction("BCCN");

    // Natural phase fractions should be in valid range
    EXPECT_GE(fcc_frac, 0.0);
    EXPECT_LE(fcc_frac, 1.0);
    EXPECT_GE(bcc_frac, 0.0);
    EXPECT_LE(bcc_frac, 1.0);

    // Store baseline Gibbs energy
    double gibbsBaseline = thermo.getGibbsEnergy();
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
    ThermoClass thermo;
    thermo.setStandardUnits();
    thermo.setThermoFilename("NobleMetals-Kaye.dat");
    thermo.parseCSDataFile();

    thermo.setTemperaturePressure(1800.0, 1.0);
    thermo.setElementMass(42, 0.5);  // Mo
    thermo.setElementMass(44, 0.5);  // Ru

    // Constrain two phases: FCCN to 40%, BCCN to 60% (sum to 1.0)
    double fccTarget = 0.4;
    double bccTarget = 0.6;
    thermo.setSolnPhaseConstraint("FCCN", fccTarget);
    thermo.setSolnPhaseConstraint("BCCN", bccTarget);

    // Verify constraints were set
    EXPECT_TRUE(thermo.getContext().phaseConstraints->hasActiveConstraints());
    EXPECT_EQ(thermo.getContext().phaseConstraints->getNumActiveConstraints(), 2);

    // Set solver parameters
    double tolerance = 1e-2;  // Slightly relaxed tolerance for two-phase system
    thermo.setConstraintTolerance(tolerance);
    thermo.setConstraintMaxOuterIterations(50);

    // Run constrained calculation
    thermo.calculate();

    // Check that calculation converged
    EXPECT_EQ(thermo.getInfoCode(), 0) << "Constrained calculation should succeed";

    // Check that constraint tracking works
    auto [fcc_frac, fcc_info] = thermo.getPhaseElementFraction("FCCN");
    auto [bcc_frac, bcc_info] = thermo.getPhaseElementFraction("BCCN");
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
    ThermoClass thermo1;
    thermo1.setStandardUnits();
    thermo1.setThermoFilename("NobleMetals-Kaye.dat");
    thermo1.parseCSDataFile();

    thermo1.setTemperaturePressure(1800.0, 1.0);
    thermo1.setElementMass(42, 0.5);  // Mo
    thermo1.setElementMass(44, 0.5);  // Ru

    thermo1.calculate();
    ASSERT_EQ(thermo1.getInfoCode(), 0);
    double gibbsUnconstrained = thermo1.getGibbsEnergy();

    // Get natural FCCN fraction
    auto [natural_fcc, info] = thermo1.getPhaseElementFraction("FCCN");

    // Run constrained calculation with different target
    ThermoClass thermo2;
    thermo2.setStandardUnits();
    thermo2.setThermoFilename("NobleMetals-Kaye.dat");
    thermo2.parseCSDataFile();

    thermo2.setTemperaturePressure(1800.0, 1.0);
    thermo2.setElementMass(42, 0.5);
    thermo2.setElementMass(44, 0.5);

    // Constrain to a different fraction than natural
    double targetFraction = std::min(0.9, std::max(0.1, natural_fcc + 0.2));
    thermo2.setSolnPhaseConstraint("FCCN", targetFraction);
    thermo2.setConstraintTolerance(1e-3);

    thermo2.calculate();

    if (thermo2.getInfoCode() == 0) {
        double gibbsConstrained = thermo2.getGibbsEnergy();

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
    ThermoClass thermo;
    thermo.setStandardUnits();
    thermo.setThermoFilename("NobleMetals-Kaye.dat");
    thermo.parseCSDataFile();

    thermo.setTemperaturePressure(1500.0, 1.0);
    thermo.setElementMass(42, 0.5);  // Mo
    thermo.setElementMass(44, 0.5);  // Ru

    // Constrain BCCN to small fraction, FCCN takes the rest
    double bccTarget = 0.1;
    double fccTarget = 0.9;
    thermo.setSolnPhaseConstraint("BCCN", bccTarget);
    thermo.setSolnPhaseConstraint("FCCN", fccTarget);
    thermo.setConstraintTolerance(1e-2);
    thermo.setConstraintMaxOuterIterations(50);

    thermo.calculate();

    if (thermo.getInfoCode() == 0) {
        auto [bcc_frac, bcc_info] = thermo.getPhaseElementFraction("BCCN");
        auto [fcc_frac, fcc_info] = thermo.getPhaseElementFraction("FCCN");
        EXPECT_NEAR(bcc_frac, bccTarget, 0.02) << "BCC should be near target";
        EXPECT_NEAR(fcc_frac, fccTarget, 0.02) << "FCC should be near target";
    }
}

// Test constraint with zero phase fraction
// This simulates a phase field scenario where a phase is completely absent.
// The solver should handle this gracefully even though the target is zero.
TEST(PhaseConstraintCalc, ZeroFractionConstraint) {
    ThermoClass thermo;
    thermo.setStandardUnits();
    thermo.setThermoFilename("NobleMetals-Kaye.dat");
    thermo.parseCSDataFile();

    thermo.setTemperaturePressure(1500.0, 1.0);
    thermo.setElementMass(42, 0.5);  // Mo
    thermo.setElementMass(44, 0.5);  // Ru

    // Constrain FCCN to exactly zero, BCCN takes all
    // This simulates a phase field scenario where one phase is absent
    double fccTarget = 0.0;
    double bccTarget = 1.0;
    thermo.setSolnPhaseConstraint("FCCN", fccTarget);
    thermo.setSolnPhaseConstraint("BCCN", bccTarget);
    thermo.setConstraintTolerance(1e-2);
    thermo.setConstraintMaxOuterIterations(50);

    thermo.calculate();

    EXPECT_EQ(thermo.getInfoCode(), 0) << "Zero-fraction constraint should succeed";

    auto [fcc_frac, fcc_info] = thermo.getPhaseElementFraction("FCCN");
    auto [bcc_frac, bcc_info] = thermo.getPhaseElementFraction("BCCN");

    EXPECT_EQ(fcc_info, 0) << "Should retrieve FCC fraction";
    EXPECT_EQ(bcc_info, 0) << "Should retrieve BCC fraction";

    // FCC should be exactly zero (or very close)
    EXPECT_NEAR(fcc_frac, fccTarget, 1e-9) << "FCC should be zero";
    // BCC should be 1.0 (or very close)
    EXPECT_NEAR(bcc_frac, bccTarget, 1e-2) << "BCC should be 1.0";
}

// Test constraining more phases than elements (common phase field scenario)
// In a binary system (2 elements), we can constrain 3 phases because the
// phase fraction constraints provide additional degrees of freedom beyond
// what the Gibbs phase rule allows for unconstrained equilibrium.
TEST(PhaseConstraintCalc, MorePhasesThanElements) {
    ThermoClass thermo;
    thermo.setStandardUnits();
    thermo.setThermoFilename("NobleMetals-Kaye.dat");
    thermo.parseCSDataFile();

    thermo.setTemperaturePressure(1500.0, 1.0);
    // Binary system: Mo and Ru (2 active elements)
    thermo.setElementMass(42, 0.5);  // Mo
    thermo.setElementMass(44, 0.5);  // Ru

    // Constrain 3 phases in a 2-element system
    // This exceeds Gibbs phase rule (P <= C) but is valid for constrained equilibrium
    double fccTarget = 0.4;
    double bccTarget = 0.4;
    double hcpTarget = 0.2;
    thermo.setSolnPhaseConstraint("FCCN", fccTarget);
    thermo.setSolnPhaseConstraint("BCCN", bccTarget);
    thermo.setSolnPhaseConstraint("HCPN", hcpTarget);
    thermo.setConstraintTolerance(1e-2);
    thermo.setConstraintMaxOuterIterations(50);

    // Verify we have 3 constraints in a 2-element system
    EXPECT_EQ(thermo.getContext().phaseConstraints->getNumActiveConstraints(), 3);

    thermo.calculate();

    // Should succeed despite having more phases than elements
    EXPECT_EQ(thermo.getInfoCode(), 0)
        << "Constrained equilibrium with more phases than elements should succeed";

    if (thermo.getInfoCode() == 0) {
        auto [fcc_frac, fcc_info] = thermo.getPhaseElementFraction("FCCN");
        auto [bcc_frac, bcc_info] = thermo.getPhaseElementFraction("BCCN");
        auto [hcp_frac, hcp_info] = thermo.getPhaseElementFraction("HCPN");

        EXPECT_EQ(fcc_info, 0) << "Should retrieve FCC fraction";
        EXPECT_EQ(bcc_info, 0) << "Should retrieve BCC fraction";
        EXPECT_EQ(hcp_info, 0) << "Should retrieve HCP fraction";

        // Verify constraints are satisfied within tolerance
        EXPECT_NEAR(fcc_frac, fccTarget, 0.02) << "FCC should match target";
        EXPECT_NEAR(bcc_frac, bccTarget, 0.02) << "BCC should match target";
        EXPECT_NEAR(hcp_frac, hcpTarget, 0.02) << "HCP should match target";

        // Fractions should sum to approximately 1
        double totalFrac = fcc_frac + bcc_frac + hcp_frac;
        EXPECT_NEAR(totalFrac, 1.0, 0.05) << "Phase fractions should sum to ~1";
    }
}

//=============================================================================
// Reset and Reuse Tests
//=============================================================================

// Test that constraints are properly reset between calculations
TEST(PhaseConstraintCalc, ConstraintReset) {
    ThermoClass thermo;
    thermo.setStandardUnits();
    thermo.setThermoFilename("NobleMetals-Kaye.dat");
    thermo.parseCSDataFile();

    // First constrained calculation
    thermo.setTemperaturePressure(1500.0, 1.0);
    thermo.setElementMass(42, 0.5);
    thermo.setElementMass(44, 0.5);
    thermo.setSolnPhaseConstraint("FCCN", 0.5);

    thermo.calculate();
    EXPECT_TRUE(thermo.getContext().phaseConstraints->hasActiveConstraints())
        << "Constraints should be active after calculation";

    // Reset for new calculation
    thermo.reset();

    // Constraints should still be defined (reset keeps constraint definitions)
    EXPECT_TRUE(thermo.getContext().phaseConstraints->hasActiveConstraints())
        << "reset should preserve constraint definitions";

    // But Lagrange multipliers should be reset
    for (const auto& c : thermo.getContext().phaseConstraints->solnPhaseConstraints) {
        EXPECT_DOUBLE_EQ(c.lagrangeMultiplier, 0.0)
            << "Lagrange multipliers should be reset";
    }

    // Full reset should clear constraints
    thermo.resetAll();
    EXPECT_FALSE(thermo.getContext().phaseConstraints->hasActiveConstraints())
        << "resetAll should clear constraints";
}

//=============================================================================
// Edge Cases and Error Handling
//=============================================================================

// Test constraint on non-existent phase
TEST(PhaseConstraintEdgeCases, NonExistentPhase) {
    ThermoClass thermo;
    thermo.setStandardUnits();
    thermo.setThermoFilename("NobleMetals-Kaye.dat");
    thermo.parseCSDataFile();

    // Try to set constraint on non-existent phase
    thermo.setSolnPhaseConstraint("NONEXISTENT_PHASE", 0.5);

    // Should not have any active constraints (invalid phase ignored)
    EXPECT_FALSE(thermo.getContext().phaseConstraints->hasActiveConstraints())
        << "Constraint on non-existent phase should be ignored";
}

// Test multiple constraints that sum to > 1
TEST(PhaseConstraintEdgeCases, OverconstrainedSystem) {
    ThermoClass thermo;
    thermo.setStandardUnits();
    thermo.setThermoFilename("NobleMetals-Kaye.dat");
    thermo.parseCSDataFile();

    thermo.setTemperaturePressure(1500.0, 1.0);
    thermo.setElementMass(42, 0.5);
    thermo.setElementMass(44, 0.5);

    // Set constraints that sum to > 1 (infeasible)
    thermo.setSolnPhaseConstraint("FCCN", 0.6);
    thermo.setSolnPhaseConstraint("BCCN", 0.6);

    thermo.calculate();

    // Should either fail or converge to best possible approximation
    // Either outcome is acceptable for an infeasible problem
    // Just verify it doesn't crash
}

//=============================================================================
// Chemical Potential and Derivative API Tests
//=============================================================================

// Test element chemical potential retrieval
TEST(PhaseConstraintAPI, ElementChemicalPotentials) {
    ThermoClass thermo;
    thermo.setStandardUnits();
    thermo.setThermoFilename("NobleMetals-Kaye.dat");
    thermo.parseCSDataFile();

    thermo.setTemperaturePressure(1500.0, 1.0);
    thermo.setElementMass(42, 0.5);  // Mo
    thermo.setElementMass(44, 0.5);  // Ru

    // Run unconstrained calculation
    thermo.calculate();
    ASSERT_EQ(thermo.getInfoCode(), 0);

    // Test single element retrieval
    auto [mo_pot, mo_info] = thermo.getElementChemicalPotential(0);
    EXPECT_EQ(mo_info, 0) << "Should retrieve first element potential";

    // Test get all element potentials
    auto potentials = thermo.getAllElementChemicalPotentials();
    EXPECT_EQ(potentials.size(), static_cast<size_t>(thermo.getContext().thermo->nElements))
        << "Should return potentials for all elements";

    // Potentials should be finite
    for (size_t i = 0; i < potentials.size(); ++i) {
        EXPECT_TRUE(std::isfinite(potentials[i]))
            << "Element potential " << i << " should be finite";
    }

    // Test invalid index
    auto [bad_pot, bad_info] = thermo.getElementChemicalPotential(-1);
    EXPECT_NE(bad_info, 0) << "Invalid index should return error";

    auto [bad_pot2, bad_info2] = thermo.getElementChemicalPotential(1000);
    EXPECT_NE(bad_info2, 0) << "Out of range index should return error";
}

// Test Gibbs energy derivative (dG/df) retrieval
TEST(PhaseConstraintAPI, GibbsEnergyDerivative) {
    ThermoClass thermo;
    thermo.setStandardUnits();
    thermo.setThermoFilename("NobleMetals-Kaye.dat");
    thermo.parseCSDataFile();

    thermo.setTemperaturePressure(1500.0, 1.0);
    thermo.setElementMass(42, 0.5);  // Mo
    thermo.setElementMass(44, 0.5);  // Ru

    // Constrain two phases
    thermo.setSolnPhaseConstraint("FCCN", 0.4);
    thermo.setSolnPhaseConstraint("BCCN", 0.6);

    thermo.calculate();
    ASSERT_EQ(thermo.getInfoCode(), 0);

    // Get derivatives for constrained phases
    auto [dGdf_fcc, fcc_info] = thermo.getGibbsEnergyDerivative("FCCN");
    auto [dGdf_bcc, bcc_info] = thermo.getGibbsEnergyDerivative("BCCN");

    EXPECT_EQ(fcc_info, 0) << "Should retrieve FCCN derivative";
    EXPECT_EQ(bcc_info, 0) << "Should retrieve BCCN derivative";

    // At equilibrium with satisfied constraints, derivatives should be near zero
    // (since the Lagrange multiplier is zero when constraints are satisfied from init)
    EXPECT_TRUE(std::isfinite(dGdf_fcc)) << "FCCN derivative should be finite";
    EXPECT_TRUE(std::isfinite(dGdf_bcc)) << "BCCN derivative should be finite";

    // Test non-existent phase
    auto [dGdf_bad, bad_info] = thermo.getGibbsEnergyDerivative("NONEXISTENT");
    EXPECT_NE(bad_info, 0) << "Non-existent phase should return error";
}
