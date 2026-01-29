# Building and Integration

[← Back to README](../README.md) | [Getting Started](getting-started.md) | [Configuration →](configuration.md)

---

How to build Thermochimica and integrate it into your project.

## Prerequisites

- **C++17 compiler**: GCC 7+, Clang 5+, or MSVC 2017+
- **CMake 3.16+**
- **Eigen3**: Auto-fetched if not found

## Building from Source

```bash
git clone https://github.com/your-repo/thermochimica.git
cd thermochimica
mkdir build && cd build
cmake ..
make -j4
```

### Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `THERMO_BUILD_TESTS` | ON | Build unit and integration tests |
| `THERMO_BUILD_EXAMPLES` | ON | Build example programs |
| `CMAKE_BUILD_TYPE` | Release | Release, Debug, or RelWithDebInfo |

**Example with options:**
```bash
cmake -DTHERMO_BUILD_TESTS=OFF -DCMAKE_BUILD_TYPE=Debug ..
```

### Running Tests

```bash
cd build
ctest --output-on-failure

# Or run specific test executable
./tests/integration_tests --gtest_filter="ThermoBasicTests.*"
```

### Running Examples

```bash
./bin/example_basic
```

---

## CMake Integration

### Option A: Subdirectory

Add Thermochimica as a subdirectory in your project:

```cmake
# CMakeLists.txt
cmake_minimum_required(VERSION 3.16)
project(MyProject)

# Add thermochimica as subdirectory
add_subdirectory(external/thermochimica)

# Your executable
add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE thermochimica)
```

### Option B: find_package (after installation)

First install Thermochimica:

```bash
cd thermochimica/build
cmake --install . --prefix /usr/local
```

Then in your project:

```cmake
cmake_minimum_required(VERSION 3.16)
project(MyProject)

find_package(Thermochimica REQUIRED)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE Thermochimica::thermochimica)
```

### Option C: FetchContent

Download automatically during CMake configure:

```cmake
cmake_minimum_required(VERSION 3.16)
project(MyProject)

include(FetchContent)
FetchContent_Declare(
    thermochimica
    GIT_REPOSITORY https://github.com/your-repo/thermochimica.git
    GIT_TAG main
)
FetchContent_MakeAvailable(thermochimica)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE thermochimica)
```

---

## Manual Linking

If not using CMake:

```bash
# Compile
g++ -std=c++17 \
    -I/path/to/thermochimica/include \
    -I/path/to/eigen3/include \
    main.cpp \
    -L/path/to/thermochimica/build \
    -lthermochimica \
    -o my_app

# Run (ensure library path is set)
export LD_LIBRARY_PATH=/path/to/thermochimica/build:$LD_LIBRARY_PATH
./my_app
```

---

## Project Structure

```
thermochimica/
├── CMakeLists.txt              # Main build file
├── include/thermochimica/      # Public headers
│   ├── Thermochimica.hpp       # Main API
│   ├── ThermoContext.hpp       # Context class
│   ├── context/                # State structures
│   ├── solver/                 # Solver headers
│   ├── parser/                 # Parser headers
│   └── util/                   # Constants, errors, tolerances
├── src/                        # Implementation
│   ├── api/                    # Public API implementation
│   ├── context/                # State management
│   ├── parser/                 # ChemSage parser
│   ├── setup/                  # Initialization
│   ├── models/                 # Thermodynamic models
│   ├── solver/                 # GEM solver
│   └── postprocess/            # Output processing
├── tests/                      # Test suites
├── examples/                   # Example programs
├── data/                       # Thermodynamic databases
└── docs/                       # Documentation
```

---

## Include Paths

After linking, use this include:

```cpp
#include <thermochimica/Thermochimica.hpp>
```

This single header includes all public API functions. For advanced usage, you can include specific headers:

```cpp
#include <thermochimica/ThermoContext.hpp>
#include <thermochimica/util/Constants.hpp>
#include <thermochimica/util/ErrorCodes.hpp>
#include <thermochimica/util/Tolerances.hpp>
```

---

## Database Files

Thermodynamic database files (`.dat`) are in the `data/` directory. You must provide the path when loading:

```cpp
// Relative path (from working directory)
Thermochimica::setThermoFilename(ctx, "data/CO.dat");

// Absolute path
Thermochimica::setThermoFilename(ctx, "/path/to/thermochimica/data/CO.dat");

// Using DATA_DIRECTORY macro (if defined in CMake)
Thermochimica::setThermoFilename(ctx, std::string(DATA_DIRECTORY) + "CO.dat");
```

---

## Compiler Flags

Recommended flags for production:

```bash
-std=c++17 -O2 -DNDEBUG
```

For debugging:

```bash
-std=c++17 -O0 -g -fsanitize=address,undefined
```

---

## Dependencies

### Eigen3

Eigen is a header-only linear algebra library. If not found on your system, CMake will fetch it automatically.

To use a specific Eigen installation:

```bash
cmake -DEigen3_DIR=/path/to/eigen3/share/eigen3/cmake ..
```

### GoogleTest (Tests Only)

Only required if building tests. Auto-fetched by CMake.

---

## Troubleshooting

### "Eigen not found"

Install Eigen or let CMake fetch it:

```bash
# Ubuntu/Debian
sudo apt install libeigen3-dev

# macOS
brew install eigen

# Or let CMake fetch (automatic)
```

### "Cannot find thermochimica.hpp"

Ensure include path is correct:

```bash
g++ -I/path/to/thermochimica/include ...
```

### "Undefined reference to Thermochimica::..."

Ensure you're linking against the library:

```bash
g++ ... -L/path/to/build -lthermochimica
```

### "Data file not found"

Use absolute path or ensure working directory is correct:

```cpp
Thermochimica::setThermoFilename(ctx, "/absolute/path/to/CO.dat");
```

---

[← Back to README](../README.md) | [Getting Started](getting-started.md) | [Configuration →](configuration.md)
