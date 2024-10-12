add_executable(sim-alnipt
    sim-alnipt.F90
)

target_link_libraries(sim-alnipt PUBLIC opencalphad octqa)

set_target_properties(sim-alnipt PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)