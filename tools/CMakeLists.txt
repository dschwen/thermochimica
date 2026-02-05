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

# Install the tools
install(TARGETS batch_validate diagnose_graphite
    RUNTIME DESTINATION bin
)
