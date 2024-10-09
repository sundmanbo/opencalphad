# CMakeLists.txt for a mixed Fortran/C++ project

cmake_minimum_required(VERSION 3.12)
project(YourProjectName LANGUAGES C CXX Fortran)

# Set the Fortran compiler
set(CMAKE_Fortran_COMPILER gfortran)
# Make sure that Fortran version is 90 or later
set(CMAKE_Fortran_STANDARD 90)

# Set C++ standard (adjust as needed)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# find openmp for C++ only
find_package(OpenMP REQUIRED COMPONENTS CXX)

# Fortran library sources
set(LIBOPENCALPHAD_SOURCES
    src/minimizer/matsmin.F90
    src/models/gtp3.F90
    src/numlib/oclablas.F90
    src/numlib/ocnum.F90
    src/numlib/minpack1.F90
    src/utilities/metlib4.F90
    examples/TQ4lib/Cpp/liboctq.F90
    examples/TQ4lib/Cpp/Scheil/liboctqisoc.F90
)

# Create libopencalphad
add_library(opencalphad ${LIBOPENCALPHAD_SOURCES})
set_target_properties(opencalphad PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Add compiler flags for internal BLAS/LAPACK implementation
target_compile_options(opencalphad PRIVATE -fno-underscoring -fno-second-underscore)

# scheil
add_executable(scheil examples/TQ4lib/Cpp/Scheil/Example_OCASI.cpp)
target_link_libraries(scheil PRIVATE 
    opencalphad
    OpenMP::OpenMP_CXX
    ${CMAKE_DL_LIBS}
)

# Ensure proper linking order
set_target_properties(scheil PROPERTIES LINK_FLAGS "-Wl,--start-group")

target_include_directories(scheil PRIVATE 
    ${CMAKE_BINARY_DIR}/examples/TQ4lib/Cpp
    ${CMAKE_BINARY_DIR}/examples/TQ4lib/Cpp/Scheil
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Ensure Fortran modules are generated before compiling C++ code
add_dependencies(scheil opencalphad)
