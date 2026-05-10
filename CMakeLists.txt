cmake_minimum_required(VERSION 3.10)
project(OpenCalphad LANGUAGES C CXX Fortran)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Set the default Fortran flags
find_package(OpenMP)

set(CMAKE_Fortran_FLAGS "-Dnotwin ${OPENMP_Fortran_FLAGS}")

# Options for building components
option(WITH_OCHELP "Build with OC Help" OFF)
option(WITH_OCPLOT "Build with OC Plot" OFF)
option(WITH_LAPACK "Use LAPACK" OFF)
option(WITH_C_BINDINGS "Export C Bindings" ON)
option(WITH_OC_EXAMPLES "Building Examples" ON)

# Define source directories
set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")
set(MODELS_DIR "${SRC_DIR}/models")
set(NUMLIB_DIR "${SRC_DIR}/numlib")
set(STEPMAP_DIR "${SRC_DIR}/stepmapplot")
set(UTILITIES_DIR "${SRC_DIR}/utilities")
set(MINIMIZER_DIR "${SRC_DIR}/minimizer")
set(USERIF_DIR "${SRC_DIR}/userif")
set(OCISOCBINDING_DIR "${CMAKE_SOURCE_DIR}/OCisoCbinding")

# Add sources for libOC library
set(LIBOC_SOURCES
    ${MODELS_DIR}/ocparam.F90
    ${UTILITIES_DIR}/metlib4.F90
    ${NUMLIB_DIR}/ocnum.F90
    ${NUMLIB_DIR}/minpack1.F90
    ${MODELS_DIR}/gtp3.F90
    ${MINIMIZER_DIR}/matsmin.F90
)

# Add conditional sources and flags
if(WITH_OCHELP)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Dlixed -Dtinyfd -Dlixhlp")
endif()

if(WITH_OCPLOT)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Dx11")
endif()

if(WITH_LAPACK)
    set(LAPACK_LIBRARIES "-llapack -lblas")
else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DNOLAPACK")
    # list(APPEND LIBOC_SOURCES ${NUMLIB_DIR}/oclablas.F90)
endif()

# Ensure correct dependency and module generation
if(NOT WITH_LAPACK)
    add_library(oc_lablas STATIC ${NUMLIB_DIR}/oclablas.F90)
    list(APPEND LIBOC_SOURCES ${NUMLIB_DIR}/ocnum.F90)
    add_library(oc_essentials STATIC ${LIBOC_SOURCES})
    target_link_libraries(oc_essentials oc_lablas ${OPENMP_Fortran_LIBRARIES})
else()
    add_library(oc_essentials STATIC ${LIBOC_SOURCES} ${NUMLIB_DIR}/ocnum.F90)
    target_link_libraries(oc_essentials ${LAPACK_LIBRARIES} ${OPENMP_Fortran_LIBRARIES})
endif()

# Add sources for the OC executable
set(OC_SOURCES
    ${STEPMAP_DIR}/smp2.F90
    ${USERIF_DIR}/pmon6.F90
    ${SRC_DIR}/pmain1.F90
)

# Create the OC executable
add_executable(OC ${OC_SOURCES})
target_link_libraries(OC oc_essentials ${LAPACK_LIBRARIES} ${OPENMP_Fortran_LIBRARIES})

if(WITH_C_BINDINGS)
    add_subdirectory(OCisoCbinding)
endif()

# Install libraries
install(TARGETS OC RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(TARGETS oc_essentials ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

if(NOT WITH_LAPACK)
    install(TARGETS oc_lablas ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

# Install extra files
install(FILES ${OCISOCBINDING_DIR}/octqc.h DESTINATION share/include)

# Install documentation
install(FILES
    ${CMAKE_SOURCE_DIR}/doc/manual/ochelp.html
    ${CMAKE_SOURCE_DIR}/doc/manual/ochelp.pdf
    DESTINATION share/doc)

if (WITH_OC_EXAMPLES)
    add_subdirectory(examples/TQ4lib)
endif()