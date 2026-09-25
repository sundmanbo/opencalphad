add_library(octqa
    liboctqa.F90
)

set_target_properties(octqa
    PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

target_link_libraries(octqa PUBLIC oc_essentials)