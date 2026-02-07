# Tools for Thermochimica validation and testing

# Batch validation tool
add_executable(batch_validate
    batch_validate.cpp
)

target_link_libraries(batch_validate
    PRIVATE
        thermochimica
)

# Graphite diagnostic tool
add_executable(diagnose_graphite
    diagnose_graphite.cpp
)

target_link_libraries(diagnose_graphite
    PRIVATE
        thermochimica
)

# Enhanced diagnostic with trace
add_executable(diagnose_with_trace
    diagnose_with_trace.cpp
)

target_link_libraries(diagnose_with_trace
    PRIVATE
        thermochimica
)

# Element order checker
add_executable(check_elements
    check_elements.cpp
)

target_link_libraries(check_elements
    PRIVATE
        thermochimica
)

# Install the tools
install(TARGETS batch_validate diagnose_graphite diagnose_with_trace
    RUNTIME DESTINATION bin
)
