# CMakeLists.txt for a mixed Fortran/C++ project

cmake_minimum_required(VERSION 3.12)
project(OpenCalphad LANGUAGES C CXX Fortran)

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
    src/models/ocparam.F90

    OCisoCbinding/liboctq.F90
    OCisoCbinding/liboctqisoc.F90
)

# Create libopencalphad
add_library(opencalphad ${LIBOPENCALPHAD_SOURCES})
set_target_properties(opencalphad PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Add compiler flags for internal BLAS/LAPACK implementation
target_compile_options(opencalphad PRIVATE -fno-underscoring -fno-second-underscore)

option(USE_EXTERNAL_LAPACK "Use external LAPACK library" OFF)

if (USE_EXTERNAL_LAPACK)
    # If LAPACK is external, do not define NOLAPACK
    target_compile_definitions(opencalphad PUBLIC -DLAPACK)
else()
    # If LAPACK is not external, define NOLAPACK
    target_compile_definitions(opencalphad PUBLIC -DNOLAPACK)
endif()

# option for including the C++ binding
option(USE_CXX_BINDING "Use C++ binding" OFF)
if (USE_CXX_BINDING)
    add_subdirectory(OCisoCbinding/cxx)
endif()

# option for including the oc6P executable
option(BUILDING_OC "Building oc6P executable" ON)
if (BUILDING_OC)
    # add executable for oc6P: src/pmain1.F90 stepmapplot/smp2.F90 stepmapplot/smp2A.F90 stepmapplot/smp2B.F90 userif/userif.F90
    add_executable(oc6p src/pmain1.F90
        src/stepmapplot/smp2.F90
        src/userif/pmon6.F90)
    target_link_libraries(oc6p PRIVATE opencalphad)
endif()
