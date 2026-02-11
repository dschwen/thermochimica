#include <gtest/gtest.h>
#include <thermochimica/ThermoClass.hpp>

using namespace Thermochimica;

// Port of test/daily/TestThermo01.F90
// Tests that an error is returned when no data file is specified

TEST(ThermoIntegration, Test01_NoDataFileSpecified) {
    ThermoClass thermo;

    // Set up conditions without specifying data file
    thermo.setStandardUnits();
    thermo.setTemperaturePressure(300.0, 1.0);
    thermo.setElementMass(6, 1.0);  // Carbon
    thermo.setElementMass(8, 1.0);  // Oxygen

    // Don't parse a data file - this should cause an error
    int result = thermo.calculate();

    // Expect error code for no data file / no species
    EXPECT_NE(result, 0);
}

// Additional integration tests would be added here
// following the same pattern as the Fortran tests

TEST(ThermoIntegration, Test_InvalidTemperature) {
    ThermoClass thermo;

    thermo.setStandardUnits();
    thermo.setTemperaturePressure(-100.0, 1.0);  // Invalid negative temperature
    thermo.setElementMass(6, 1.0);

    // Would need database loaded first
    // int result = thermo.calculate();
    // EXPECT_EQ(result, ErrorCode::kTemperatureOutOfRange);
}

TEST(ThermoIntegration, Test_InvalidPressure) {
    ThermoClass thermo;

    thermo.setStandardUnits();
    thermo.setTemperaturePressure(1000.0, -1.0);  // Invalid negative pressure
    thermo.setElementMass(6, 1.0);

    // Would need database loaded first
    // int result = thermo.calculate();
    // EXPECT_EQ(result, ErrorCode::kPressureOutOfRange);
}

TEST(ThermoIntegration, Test_NoComposition) {
    ThermoClass thermo;

    thermo.setStandardUnits();
    thermo.setTemperaturePressure(1000.0, 1.0);
    // Don't set any element masses

    // Would need database loaded first
    // int result = thermo.calculate();
    // EXPECT_EQ(result, ErrorCode::kCompositionOutOfRange);
}
