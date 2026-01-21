add_executable(tq1-crfe 
    TQ1-crfe.F90
)

target_link_libraries(tq1-crfe PRIVATE oc_essentials octqa)