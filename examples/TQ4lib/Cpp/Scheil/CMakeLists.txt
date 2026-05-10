# scheil
add_executable(scheil 
    Example_OCASI.cpp 
)

target_link_libraries(scheil PRIVATE 
    oc_essentials
    OPENCALPHAD
    OpenMP::OpenMP_CXX
    ${CMAKE_DL_LIBS}
)

# Ensure proper linking order
set_target_properties(scheil PROPERTIES LINK_FLAGS "-Wl,--start-group")

target_include_directories(scheil PRIVATE 
    ${CMAKE_BINARY_DIR}/
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Ensure Fortran modules are generated before compiling C++ code
add_dependencies(scheil oc_essentials)
