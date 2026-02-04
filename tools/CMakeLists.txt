# Tools for Thermochimica validation and testing

# Batch validation tool
add_executable(batch_validate
    batch_validate.cpp
)

target_link_libraries(batch_validate
    PRIVATE
        thermochimica
)

# Install the tool
install(TARGETS batch_validate
    RUNTIME DESTINATION bin
)
