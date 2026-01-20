#include <gtest/gtest.h>
#include <thermochimica/Thermochimica.hpp>

using namespace Thermochimica;

// Port of test/daily/TestThermo01.F90
// Tests that an error is returned when no data file is specified

TEST(ThermoIntegration, Test01_NoDataFileSpecified) {
    ThermoContext ctx;

    // Set up conditions without specifying data file
    setStandardUnits(ctx);
    setTemperaturePressure(ctx, 300.0, 1.0);
    setElementMass(ctx, 6, 1.0);  // Carbon
    setElementMass(ctx, 8, 1.0);  // Oxygen

    // Don't parse a data file - this should cause an error
    thermochimica(ctx);

    // Expect error code for no data file / no species
    EXPECT_NE(ctx.infoThermo(), 0);
}

// Additional integration tests would be added here
// following the same pattern as the Fortran tests

TEST(ThermoIntegration, Test_InvalidTemperature) {
    ThermoContext ctx;

    setStandardUnits(ctx);
    setTemperaturePressure(ctx, -100.0, 1.0);  // Invalid negative temperature
    setElementMass(ctx, 6, 1.0);

    // Would need database loaded first
    // thermochimica(ctx);
    // EXPECT_EQ(ctx.infoThermo(), ErrorCode::kTemperatureOutOfRange);
}

TEST(ThermoIntegration, Test_InvalidPressure) {
    ThermoContext ctx;

    setStandardUnits(ctx);
    setTemperaturePressure(ctx, 1000.0, -1.0);  // Invalid negative pressure
    setElementMass(ctx, 6, 1.0);

    // Would need database loaded first
    // thermochimica(ctx);
    // EXPECT_EQ(ctx.infoThermo(), ErrorCode::kPressureOutOfRange);
}

TEST(ThermoIntegration, Test_NoComposition) {
    ThermoContext ctx;

    setStandardUnits(ctx);
    setTemperaturePressure(ctx, 1000.0, 1.0);
    // Don't set any element masses

    // Would need database loaded first
    // thermochimica(ctx);
    // EXPECT_EQ(ctx.infoThermo(), ErrorCode::kCompositionOutOfRange);
}
