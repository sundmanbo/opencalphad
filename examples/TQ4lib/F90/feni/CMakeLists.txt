add_executable(tq2-feni
    TQ2-feni.F90
)

target_link_libraries(tq2-feni PUBLIC oc_essentials octqa)

# set_target_properties(tq2-feni PROPERTIES
#     Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
# )