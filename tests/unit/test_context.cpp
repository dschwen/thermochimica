#include <gtest/gtest.h>
#include <thermochimica/ThermoContext.hpp>
#include <thermochimica/Thermochimica.hpp>
#include <thermochimica/ThermoClass.hpp>
#include <thermochimica/util/ErrorCodes.hpp>

using namespace Thermochimica;

TEST(ContextTest, Creation) {
    ThermoContext ctx;

    EXPECT_NE(ctx.thermo, nullptr);
    EXPECT_NE(ctx.io, nullptr);
    EXPECT_NE(ctx.gem, nullptr);
    EXPECT_NE(ctx.parser, nullptr);
    EXPECT_NE(ctx.submin, nullptr);
    EXPECT_NE(ctx.ctz, nullptr);
    EXPECT_NE(ctx.reinit, nullptr);
}

TEST(ContextTest, InitialState) {
    ThermoContext ctx;

    EXPECT_EQ(ctx.infoThermo(), 0);
    EXPECT_TRUE(ctx.isSuccess());
    EXPECT_FALSE(ctx.isDatabaseLoaded());
}

TEST(ContextTest, SetUnits) {
    ThermoClass thermo;
    thermo.setStandardUnits();

    auto& ctx = thermo.getContext();
    EXPECT_EQ(ctx.io->cInputUnitTemperature, "K");
    EXPECT_EQ(ctx.io->cInputUnitPressure, "atm");
    EXPECT_EQ(ctx.io->cInputUnitMass, "moles");
}

TEST(ContextTest, SetTemperaturePressure) {
    ThermoClass thermo;
    thermo.setTemperaturePressure(1000.0, 1.0);

    auto& ctx = thermo.getContext();
    // Values are stored as raw input until initialize() converts them
    EXPECT_DOUBLE_EQ(ctx.io->dTemperatureInput, 1000.0);
    EXPECT_DOUBLE_EQ(ctx.io->dPressureInput, 1.0);
    EXPECT_FALSE(ctx.io->bTemperatureConverted);
    EXPECT_FALSE(ctx.io->bPressureConverted);
}

TEST(ContextTest, SetElementMass) {
    ThermoClass thermo;
    thermo.setElementMass(6, 1.0);   // Carbon
    thermo.setElementMass(8, 2.0);   // Oxygen

    auto& ctx = thermo.getContext();
    EXPECT_DOUBLE_EQ(ctx.io->dElementMass[6], 1.0);
    EXPECT_DOUBLE_EQ(ctx.io->dElementMass[8], 2.0);
}

TEST(ContextTest, Reset) {
    ThermoClass thermo;
    thermo.setTemperaturePressure(1000.0, 1.0);
    thermo.setElementMass(6, 1.0);

    thermo.reset();

    auto& ctx = thermo.getContext();
    EXPECT_EQ(ctx.infoThermo(), 0);
    EXPECT_DOUBLE_EQ(ctx.io->dElementMass[6], 0.0);
}

TEST(TolerancesTest, Defaults) {
    Tolerances tol;

    EXPECT_GT(tol[kTolMassBalance], 0.0);
    EXPECT_GT(tol[kTolChemPotential], 0.0);
    EXPECT_GT(tol[kTolConvergence], 0.0);
}

TEST(ErrorCodeTest, Messages) {
    EXPECT_STREQ(ErrorCode::getMessage(0), "Success");
    EXPECT_STREQ(ErrorCode::getMessage(ErrorCode::kNoDataFileSpecified),
                 "No thermodynamic data file specified");
}

TEST(AtomicNumberTest, Lookup) {
    EXPECT_EQ(getAtomicNumber("H"), 1);
    EXPECT_EQ(getAtomicNumber("C"), 6);
    EXPECT_EQ(getAtomicNumber("O"), 8);
    EXPECT_EQ(getAtomicNumber("Fe"), 26);
    EXPECT_EQ(getAtomicNumber("U"), 92);
    EXPECT_EQ(getAtomicNumber("Unknown"), -1);
}

TEST(ElementSymbolTest, Lookup) {
    EXPECT_EQ(getElementSymbol(1), "H");
    EXPECT_EQ(getElementSymbol(6), "C");
    EXPECT_EQ(getElementSymbol(8), "O");
    EXPECT_EQ(getElementSymbol(26), "Fe");
    EXPECT_EQ(getElementSymbol(92), "U");
    EXPECT_EQ(getElementSymbol(0), "");
    EXPECT_EQ(getElementSymbol(200), "");
}
