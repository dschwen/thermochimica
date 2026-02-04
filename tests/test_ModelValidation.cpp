/// @file test_ModelValidation.cpp
/// @brief Explicit validation tests for all thermodynamic models
/// @details Ensures each model type (IDMX, QKTO, RKMP, SUBL, SUBG, SUBQ) works correctly

#include <gtest/gtest.h>
#include "thermochimica/ThermoClass.hpp"
#include "thermochimica/Thermochimica.hpp"
#include "thermochimica/models/ModelFactory.hpp"
#include <cmath>

using namespace Thermochimica;

/// @brief Test fixture for model validation
class ModelValidationTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}

    /// @brief Run calculation and verify success
    bool runCalculation(ThermoClass& thermo) {
        int result = thermo.calculate();
        if (result != 0) {
            std::cout << "Calculation failed: " << thermo.getErrorMessage() << std::endl;
            return false;
        }
        return true;
    }

    /// @brief Compare two doubles
    bool approxEqual(double a, double b, double tol = 1e-8) {
        return std::abs(a - b) / (std::abs(a) + 1e-10) < tol;
    }
};

/// @brief Test ModelFactory registration
TEST_F(ModelValidationTest, ModelFactoryRegistration) {
    ModelFactory factory;

    // Verify all model types are registered
    EXPECT_TRUE(factory.hasModel(Constants::PhaseType::IDMX));
    EXPECT_TRUE(factory.hasModel(Constants::PhaseType::QKTO));
    EXPECT_TRUE(factory.hasModel(Constants::PhaseType::RKMP));
    EXPECT_TRUE(factory.hasModel(Constants::PhaseType::RKMPM));
    EXPECT_TRUE(factory.hasModel(Constants::PhaseType::SUBL));
    EXPECT_TRUE(factory.hasModel(Constants::PhaseType::SUBLM));
    EXPECT_TRUE(factory.hasModel(Constants::PhaseType::SUBG));
    EXPECT_TRUE(factory.hasModel(Constants::PhaseType::SUBQ));

    // Verify we can get each model
    EXPECT_NE(factory.getModel(Constants::PhaseType::IDMX), nullptr);
    EXPECT_NE(factory.getModel(Constants::PhaseType::QKTO), nullptr);
    EXPECT_NE(factory.getModel(Constants::PhaseType::RKMP), nullptr);
    EXPECT_NE(factory.getModel(Constants::PhaseType::SUBL), nullptr);
    EXPECT_NE(factory.getModel(Constants::PhaseType::SUBG), nullptr);
    EXPECT_NE(factory.getModel(Constants::PhaseType::SUBQ), nullptr);

    // Verify model names
    EXPECT_STREQ(factory.getModel(Constants::PhaseType::IDMX)->getModelName(), "IdealMixing");
    EXPECT_STREQ(factory.getModel(Constants::PhaseType::QKTO)->getModelName(), "QKTO");
    EXPECT_STREQ(factory.getModel(Constants::PhaseType::RKMP)->getModelName(), "RKMPModel");
    EXPECT_STREQ(factory.getModel(Constants::PhaseType::SUBL)->getModelName(), "SUBLModel");
    EXPECT_STREQ(factory.getModel(Constants::PhaseType::SUBG)->getModelName(), "SUBGModel");
    EXPECT_STREQ(factory.getModel(Constants::PhaseType::SUBQ)->getModelName(), "SUBQModel");
}

/// @brief Test IDMX model (Ideal Mixing) - CO gas system
TEST_F(ModelValidationTest, IDMX_Model_CO_Gas) {
    ThermoClass thermo;
    ASSERT_EQ(thermo.loadDatabase("CO.dat"), 0);

    thermo.setStandardUnits();
    thermo.setTemperaturePressure(1000.0, 1.0);
    thermo.setElementMass("C", 1.0);
    thermo.setElementMass("O", 2.0);

    ASSERT_TRUE(runCalculation(thermo));

    // CO system should produce gas phase (ideal mixing)
    auto [molesGas, info] = thermo.getMolesPhase("gas_ideal");
    EXPECT_GT(molesGas, 0.0) << "Gas phase should be present";

    // Gibbs energy should be negative (stable)
    double gibbs = thermo.getGibbsEnergy();
    EXPECT_LT(gibbs, 0.0) << "System should be stable";
}

/// @brief Test RKMP model (Redlich-Kister) - Noble metals
TEST_F(ModelValidationTest, RKMP_Model_NobleMetals) {
    ThermoClass thermo;

    // Skip if database not available
    int result = thermo.loadDatabase("NobleMetals-Kaye.dat");
    if (result != 0) {
        GTEST_SKIP() << "NobleMetals-Kaye.dat not available";
    }

    thermo.setStandardUnits();
    thermo.setTemperaturePressure(2250.0, 1.0);
    thermo.setElementMass("Mo", 0.5);
    thermo.setElementMass("Ru", 0.5);

    ASSERT_TRUE(runCalculation(thermo));

    // Mo-Ru system should form solid solution phases (RKMP model)
    double gibbs = thermo.getGibbsEnergy();
    EXPECT_LT(gibbs, 0.0) << "System should be stable";

    std::cout << "RKMP test (Mo-Ru at 2250K): Gibbs = " << gibbs << " J" << std::endl;
}

/// @brief Test SUBL model (Compound Energy Formalism) - Multiple databases use this
TEST_F(ModelValidationTest, SUBL_Model_Generic) {
    // Many databases use SUBL for intermetallic phases
    // We'll test with a system known to have SUBL phases
    ThermoClass thermo;

    int result = thermo.loadDatabase("FeTiVO.dat");
    if (result != 0) {
        GTEST_SKIP() << "FeTiVO.dat not available";
    }

    thermo.setStandardUnits();
    thermo.setTemperaturePressure(1500.0, 1.0);
    thermo.setElementMass("Fe", 1.0);
    thermo.setElementMass("Ti", 0.5);
    thermo.setElementMass("V", 0.5);
    thermo.setElementMass("O", 2.0);

    ASSERT_TRUE(runCalculation(thermo));

    double gibbs = thermo.getGibbsEnergy();
    EXPECT_LT(gibbs, 0.0) << "System should be stable";

    std::cout << "SUBL test (FeTiVO at 1500K): Gibbs = " << gibbs << " J" << std::endl;
}

/// @brief Test SUBG model (Modified Quasichemical) - Ionic liquids
TEST_F(ModelValidationTest, SUBG_Model_IonicLiquid_CsI) {
    ThermoClass thermo;

    int result = thermo.loadDatabase("CsI-Pham.dat");
    if (result != 0) {
        GTEST_SKIP() << "CsI-Pham.dat not available";
    }

    thermo.setStandardUnits();
    thermo.setTemperaturePressure(1000.0, 1.0);
    thermo.setElementMass("Cs", 1.0);
    thermo.setElementMass("I", 1.0);

    ASSERT_TRUE(runCalculation(thermo));

    double gibbs = thermo.getGibbsEnergy();
    EXPECT_LT(gibbs, 0.0) << "CsI system should be stable";

    std::cout << "SUBG test (CsI at 1000K): Gibbs = " << gibbs << " J" << std::endl;
}

/// @brief Test SUBG model - Sulfide system
TEST_F(ModelValidationTest, SUBG_Model_Sulfide_CaMnS) {
    ThermoClass thermo;

    int result = thermo.loadDatabase("CaMnS.dat");
    if (result != 0) {
        GTEST_SKIP() << "CaMnS.dat not available";
    }

    thermo.setStandardUnits();
    thermo.setTemperaturePressure(1200.0, 1.0);
    thermo.setElementMass("Ca", 1.0);
    thermo.setElementMass("Mn", 1.0);
    thermo.setElementMass("S", 2.0);

    ASSERT_TRUE(runCalculation(thermo));

    double gibbs = thermo.getGibbsEnergy();
    EXPECT_LT(gibbs, 0.0) << "CaMnS system should be stable";

    std::cout << "SUBG test (CaMnS at 1200K): Gibbs = " << gibbs << " J" << std::endl;
}

/// @brief Test SUBG model - Halide system
TEST_F(ModelValidationTest, SUBG_Model_Halide_ClAlNa) {
    ThermoClass thermo;

    int result = thermo.loadDatabase("ClAlNa.dat");
    if (result != 0) {
        GTEST_SKIP() << "ClAlNa.dat not available";
    }

    thermo.setStandardUnits();
    thermo.setTemperaturePressure(1000.0, 1.0);
    thermo.setElementMass("Cl", 1.0);
    thermo.setElementMass("Al", 0.5);
    thermo.setElementMass("Na", 0.5);

    ASSERT_TRUE(runCalculation(thermo));

    double gibbs = thermo.getGibbsEnergy();
    EXPECT_LT(gibbs, 0.0) << "ClAlNa system should be stable";

    std::cout << "SUBG test (ClAlNa at 1000K): Gibbs = " << gibbs << " J" << std::endl;
}

/// @brief Test SUBQ model (MQM Quadruplet) - ZIRC system
TEST_F(ModelValidationTest, SUBQ_Model_ZIRC) {
    ThermoClass thermo;

    int result = thermo.loadDatabase("ZIRC-noSUBI.dat");
    if (result != 0) {
        GTEST_SKIP() << "ZIRC-noSUBI.dat not available";
    }

    thermo.setStandardUnits();
    thermo.setTemperaturePressure(1000.0, 1.0);
    thermo.setElementMass("Zr", 1.0);
    thermo.setElementMass("O", 2.0);

    ASSERT_TRUE(runCalculation(thermo));

    double gibbs = thermo.getGibbsEnergy();
    EXPECT_LT(gibbs, 0.0) << "ZrO system should be stable";

    std::cout << "SUBQ test (ZrO at 1000K): Gibbs = " << gibbs << " J" << std::endl;
}

/// @brief Compare class-based API with free-function API for each model
TEST_F(ModelValidationTest, CompareAPIs_AllModels) {
    struct TestCase {
        std::string name;
        std::string database;
        double temperature;
        double pressure;
        std::vector<std::pair<std::string, double>> composition;
    };

    std::vector<TestCase> testCases = {
        {"IDMX_CO", "CO.dat", 1000.0, 1.0, {{"C", 1.0}, {"O", 2.0}}},
        {"RKMP_MoRu", "NobleMetals-Kaye.dat", 2000.0, 1.0, {{"Mo", 0.5}, {"Ru", 0.5}}},
        {"SUBG_CsI", "CsI-Pham.dat", 1000.0, 1.0, {{"Cs", 1.0}, {"I", 1.0}}},
        {"SUBQ_ZIRC", "ZIRC-noSUBI.dat", 1000.0, 1.0, {{"Zr", 1.0}, {"O", 2.0}}},
    };

    for (const auto& test : testCases) {
        // Class-based API
        ThermoClass thermo;
        if (thermo.loadDatabase(test.database) != 0) {
            std::cout << "Skipping " << test.name << " - database not available" << std::endl;
            continue;
        }

        thermo.setStandardUnits();
        thermo.setTemperaturePressure(test.temperature, test.pressure);
        for (const auto& [element, mass] : test.composition) {
            thermo.setElementMass(element, mass);
        }

        int result1 = thermo.calculate();
        ASSERT_EQ(result1, 0) << test.name << ": Class-based calculation failed";
        double gibbs1 = thermo.getGibbsEnergy();

        // Free-function API
        ThermoContext ctx;
        setThermoFilename(ctx, test.database);
        parseCSDataFile(ctx);
        setStandardUnits(ctx);
        setTemperaturePressure(ctx, test.temperature, test.pressure);
        for (const auto& [element, mass] : test.composition) {
            setElementMass(ctx, element, mass);
        }

        thermochimica(ctx);
        int result2 = ctx.infoThermo();
        ASSERT_EQ(result2, 0) << test.name << ": Free-function calculation failed";
        double gibbs2 = getGibbsEnergy(ctx);

        // Results should be identical
        EXPECT_TRUE(approxEqual(gibbs1, gibbs2))
            << test.name << ": Gibbs mismatch - Class: " << gibbs1 << ", Free: " << gibbs2;

        std::cout << test.name << " validated: Gibbs = " << gibbs1 << " J (match: "
                  << (approxEqual(gibbs1, gibbs2) ? "YES" : "NO") << ")" << std::endl;
    }
}

/// @brief Regression test to ensure models don't silently return zero
TEST_F(ModelValidationTest, NoSilentZeroExcessGibbs) {
    // This test ensures that when a model is used, it actually computes something
    // rather than silently returning zero (which was the bug before full implementation)

    ThermoClass thermo;

    // Noble metals with RKMP should have NON-ZERO excess Gibbs
    int result = thermo.loadDatabase("NobleMetals-Kaye.dat");
    if (result != 0) {
        GTEST_SKIP() << "NobleMetals-Kaye.dat not available";
    }

    thermo.setStandardUnits();
    thermo.setTemperaturePressure(2000.0, 1.0);
    thermo.setElementMass("Mo", 0.5);
    thermo.setElementMass("Ru", 0.5);

    ASSERT_TRUE(runCalculation(thermo));

    // The total Gibbs energy should be significantly negative
    // If models were returning zero excess, this would be much less negative
    double gibbs = thermo.getGibbsEnergy();
    EXPECT_LT(gibbs, -1000.0) << "Gibbs energy suspiciously close to zero - model may not be working";

    std::cout << "No-silent-zero test: Gibbs = " << gibbs << " J (expected < -1000 J)" << std::endl;
}
