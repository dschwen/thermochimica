cmake_minimum_required(VERSION 3.16)
project(Thermochimica VERSION 2.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Options
option(THERMO_BUILD_TESTS "Build tests" ON)
option(THERMO_BUILD_EXAMPLES "Build examples" ON)
option(THERMO_BUILD_TOOLS "Build validation tools" ON)

# Find Eigen3 (or fetch it if not found)
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    message(STATUS "Eigen3 not found, fetching from git...")
    include(FetchContent)
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(eigen)
endif()

# Define data directory - use STRING to preserve trailing slash
set(THERMO_DATA_DIRECTORY "${CMAKE_SOURCE_DIR}/data/" CACHE STRING "Path to thermodynamic data files")
add_compile_definitions(DATA_DIRECTORY="${THERMO_DATA_DIRECTORY}")

# Compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Main library sources
set(THERMO_SOURCES
    # Context
    src/context/ThermoContext.cpp
    src/context/ThermoState.cpp
    src/context/ThermoIO.cpp
    src/context/GEMState.cpp
    src/context/ParserState.cpp
    src/context/SubMinState.cpp
    src/context/CTZState.cpp
    src/context/ReinitState.cpp
    src/context/PhaseConstraints.cpp

    # Parser
    src/parser/ChemSageParser.cpp

    # Setup
    src/setup/InitThermo.cpp
    src/setup/CheckThermoInput.cpp

    # Models
    src/models/ExcessGibbs.cpp
    src/models/IdealMixing.cpp
    src/models/QKTO.cpp
    src/models/RKMP.cpp
    src/models/SUBL.cpp
    src/models/SUBG.cpp
    src/models/SUBQ.cpp
    # Model classes (OO architecture)
    src/models/IdealMixingModel.cpp
    src/models/QKTOModel.cpp
    src/models/RKMPModel.cpp
    src/models/SUBLModel.cpp
    src/models/SUBGModel.cpp
    src/models/SUBQModel.cpp
    src/models/ModelFactory.cpp

    # Solver
    src/solver/GEMSolver.cpp
    src/solver/StandardGEMSolver.cpp
    src/solver/ConstrainedGEMSolver.cpp
    src/solver/GEMNewton.cpp
    src/solver/NewtonSolver.cpp
    src/solver/GEMLineSearch.cpp
    src/solver/WolfeLineSearch.cpp
    src/solver/LevelingSolver.cpp
    src/solver/Subminimization.cpp
    src/solver/PhaseAssemblage.cpp
    src/solver/PhaseAssemblageManager.cpp
    src/solver/CheckConvergence.cpp
    src/solver/ConstrainedGEM.cpp

    # Post-process
    src/postprocess/PostProcess.cpp
    src/postprocess/PrintResults.cpp
    src/postprocess/HeatCapacity.cpp

    # API
    src/api/Thermochimica.cpp
    src/api/ThermoClass.cpp
    src/api/CouplingAPI.cpp
)

# Main library
add_library(thermochimica ${THERMO_SOURCES})

target_include_directories(thermochimica
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(thermochimica PUBLIC Eigen3::Eigen)

# Alias for cleaner usage
add_library(Thermochimica::thermochimica ALIAS thermochimica)

# Tests
if(THERMO_BUILD_TESTS)
    enable_testing()
    find_package(GTest QUIET)
    if(NOT GTest_FOUND)
        message(STATUS "GTest not found, fetching from git...")
        include(FetchContent)
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG v1.14.0
            GIT_SHALLOW TRUE
        )
        # For Windows: Prevent overriding the parent project's compiler/linker settings
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(googletest)
    endif()
    add_subdirectory(tests)
endif()

# Examples
if(THERMO_BUILD_EXAMPLES)
    add_executable(example_basic examples/basic_usage.cpp)
    target_link_libraries(example_basic PRIVATE thermochimica)
endif()

# Tools
if(THERMO_BUILD_TOOLS)
    add_subdirectory(tools)
endif()

# Installation
include(GNUInstallDirs)

install(TARGETS thermochimica
    EXPORT ThermochimicaTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT ThermochimicaTargets
    FILE ThermochimicaTargets.cmake
    NAMESPACE Thermochimica::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Thermochimica
)

# Package config
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_SOURCE_DIR}/cmake/ThermochimicaConfig.cmake.in
    ${CMAKE_BINARY_DIR}/ThermochimicaConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Thermochimica
)

write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/ThermochimicaConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_BINARY_DIR}/ThermochimicaConfig.cmake
    ${CMAKE_BINARY_DIR}/ThermochimicaConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Thermochimica
)
