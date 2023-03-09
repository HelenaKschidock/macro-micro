#copy all source code etc to build directory
add_custom_target(copy_micro_heat_all ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(pybind11)
pybind11_add_module(micro_sim micro_sim.cpp)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
#target_compile_definitions(micro_sim
#                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
target_link_libraries(micro_sim PRIVATE dumux-precice)

# add a symlink for each input file
add_input_file_links()
