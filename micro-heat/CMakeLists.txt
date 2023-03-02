#copy all source code etc to build directory
add_custom_target(copy_micro_heat_all ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_micro_heat EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_micro_heat PUBLIC "ENABLEMONOLITHIC=0")
target_link_libraries(test_micro_heat PRIVATE dumux-precice)

add_executable(run_micro_heat_sims EXCLUDE_FROM_ALL micro_sim.cpp)
target_compile_definitions(run_micro_heat_sims PUBLIC "ENABLEMONOLITHIC=0")
target_link_libraries(run_micro_heat_sims PRIVATE dumux-precice)

# add a symlink for each input file
add_input_file_links()
