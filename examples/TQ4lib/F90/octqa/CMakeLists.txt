add_library(octqa
    liboctqa.F90
)

set_target_properties(octqa PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_link_libraries(octqa PUBLIC opencalphad)