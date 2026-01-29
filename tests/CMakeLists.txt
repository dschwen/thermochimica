# Test configuration

include(GoogleTest)

# Handle both find_package and FetchContent GTest
if(TARGET GTest::gtest)
    set(GTEST_LIBS GTest::gtest GTest::gtest_main)
else()
    set(GTEST_LIBS gtest gtest_main)
endif()

# Unit tests
add_executable(unit_tests
    unit/test_context.cpp
)

target_link_libraries(unit_tests
    PRIVATE
        thermochimica
        ${GTEST_LIBS}
)

gtest_discover_tests(unit_tests)

# Integration tests (ported from Fortran)
add_executable(integration_tests
    integration/TestThermo01.cpp
    integration/TestThermo.cpp
)

target_link_libraries(integration_tests
    PRIVATE
        thermochimica
        ${GTEST_LIBS}
)

gtest_discover_tests(integration_tests)
