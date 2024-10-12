add_executable(sim-alnipt
    sim-alnipt.F90
)

target_link_libraries(sim-alnipt PUBLIC oc_essentials octqa)
