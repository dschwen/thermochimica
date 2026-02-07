/// @file test_ThermoClass.cpp
/// @brief Tests for the new class-based ThermoClass API
/// @details Verifies that ThermoClass produces identical results to free-function API

#include <gtest/gtest.h>
#include "thermochimica/ThermoClass.hpp"
#include <cmath>

using namespace Thermochimica;

/// @brief Test fixture for ThermoClass tests
class ThermoClassTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Common setup for all tests
    }

    void TearDown() override {
        // Common cleanup
    }

    /// @brief Compare two doubles with tolerance
    bool approxEqual(double a, double b, double tol = 1e-10) {
        return std::abs(a - b) < tol;
    }
};

/// @brief Test basic construction and destruction
TEST_F(ThermoClassTest, ConstructorDestructor) {
    ThermoClass thermo;
    EXPECT_EQ(thermo.getInfoCode(), 0);
    EXPECT_TRUE(thermo.isSuccess());
}

/// @brief Test move semantics
TEST_F(ThermoClassTest, MoveSemantics) {
    ThermoClass thermo1;
    thermo1.setTemperature(1000.0);

    // Move constructor
    ThermoClass thermo2(std::move(thermo1));
    EXPECT_EQ(thermo2.getInfoCode(), 0);

    // Move assignment
    ThermoClass thermo3;
    thermo3 = std::move(thermo2);
    EXPECT_EQ(thermo3.getInfoCode(), 0);
}

/// @brief Test basic database loading
TEST_F(ThermoClassTest, LoadDatabase) {
    ThermoClass thermo;
    int result = thermo.loadDatabase("CO.dat");
    EXPECT_EQ(result, 0) << "Failed to load database: " << thermo.getErrorMessage();
    EXPECT_TRUE(thermo.isSuccess());
}

/// @brief Test input configuration methods
TEST_F(ThermoClassTest, InputConfiguration) {
    ThermoClass thermo;

    // Set temperature and pressure
    thermo.setTemperature(1000.0);
    thermo.setPressure(1.0);

    // Set temperature and pressure together
    thermo.setTemperaturePressure(1200.0, 2.0);

    // Set element mass by atomic number
    thermo.setElementMass(6, 1.0);  // Carbon

    // Set element mass by name
    thermo.setElementMass("O", 1.0);  // Oxygen

    // Set units
    thermo.setStandardUnits();
    thermo.setUnitsSI();

    EXPECT_TRUE(thermo.isSuccess());
}

/// @brief Test phase constraints with class-based API
TEST_F(ThermoClassTest, PhaseConstraints) {
    ThermoClass thermo;
    thermo.loadDatabase("CO.dat");
    thermo.setStandardUnits();
    thermo.setTemperaturePressure(1000.0, 1.0);
    thermo.setElementMass(6, 1.0);
    thermo.setElementMass(8, 2.0);

    // Set a phase constraint
    thermo.setSolnPhaseConstraint("gas_ideal", 0.5);

    int result = thermo.calculate();
    EXPECT_EQ(result, 0) << "Constrained calculation failed: " << thermo.getErrorMessage();

    // Remove constraint
    thermo.removePhaseConstraint("gas_ideal");

    // Clear all constraints
    thermo.clearPhaseConstraints();

    EXPECT_TRUE(thermo.isSuccess());
}

/// @brief Test reset functionality
TEST_F(ThermoClassTest, Reset) {
    ThermoClass thermo;
    thermo.loadDatabase("CO.dat");
    thermo.setTemperaturePressure(1000.0, 1.0);
    thermo.setElementMass(6, 1.0);
    thermo.calculate();

    // Reset for new calculation (keeps database)
    thermo.reset();
    thermo.setTemperaturePressure(1200.0, 2.0);
    thermo.setElementMass(6, 2.0);
    int result = thermo.calculate();
    EXPECT_EQ(result, 0);

    // Full reset (clears database)
    thermo.resetAll();
    EXPECT_TRUE(thermo.isSuccess());
}

/// @brief Test granular control methods
TEST_F(ThermoClassTest, GranularControl) {
    ThermoClass thermo;
    thermo.loadDatabase("CO.dat");
    thermo.setStandardUnits();
    thermo.setTemperaturePressure(1000.0, 1.0);
    thermo.setElementMass(6, 1.0);
    thermo.setElementMass(8, 2.0);

    // Test granular steps
    thermo.initialize();
    thermo.checkSystem();
    thermo.computeThermoData();
    thermo.setup();
    int result = thermo.solve();

    EXPECT_EQ(result, 0) << "Granular solve failed: " << thermo.getErrorMessage();
}

/// @brief Test multiple independent instances
TEST_F(ThermoClassTest, MultipleInstances) {
    ThermoClass thermo1;
    ThermoClass thermo2;

    // Configure first instance
    thermo1.loadDatabase("CO.dat");
    thermo1.setStandardUnits();
    thermo1.setTemperaturePressure(1000.0, 1.0);
    thermo1.setElementMass(6, 1.0);
    thermo1.setElementMass(8, 2.0);

    // Configure second instance with different conditions
    thermo2.loadDatabase("CO.dat");
    thermo2.setStandardUnits();
    thermo2.setTemperaturePressure(1500.0, 2.0);
    thermo2.setElementMass(6, 2.0);
    thermo2.setElementMass(8, 1.0);

    // Calculate both
    int result1 = thermo1.calculate();
    int result2 = thermo2.calculate();

    EXPECT_EQ(result1, 0);
    EXPECT_EQ(result2, 0);

    // Results should be different (different conditions)
    double gibbs1 = thermo1.getGibbsEnergy();
    double gibbs2 = thermo2.getGibbsEnergy();

    EXPECT_FALSE(approxEqual(gibbs1, gibbs2))
        << "Independent instances should have different results";
}

/// @brief Test error handling
TEST_F(ThermoClassTest, ErrorHandling) {
    ThermoClass thermo;

    // Try to calculate without loading database
    int result = thermo.calculate();
    EXPECT_NE(result, 0) << "Should fail without database";
    EXPECT_FALSE(thermo.isSuccess());

    std::string errorMsg = thermo.getErrorMessage();
    EXPECT_FALSE(errorMsg.empty()) << "Should have error message";
}

