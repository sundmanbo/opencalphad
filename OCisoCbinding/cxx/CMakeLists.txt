# scheil
add_executable(test_binding 
    Example_OCASI.cpp 
)

target_link_libraries(test_binding PRIVATE 
    opencalphad
    oc_binding
    OpenMP::OpenMP_CXX
    ${CMAKE_DL_LIBS}
)

# Ensure proper linking order
set_target_properties(test_binding PROPERTIES LINK_FLAGS "-Wl,--start-group")

target_include_directories(test_binding PRIVATE 
    ${CMAKE_BINARY_DIR}/OCisoCbinding/cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Ensure Fortran modules are generated before compiling C++ code
add_dependencies(test_binding opencalphad)
