#copy all source code etc to build directory
add_custom_target(copy_micro_heat_all ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_executable(test_micro_heat EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_micro_heat PUBLIC "ENABLEMONOLITHIC=0")
target_link_libraries(test_micro_heat PRIVATE dumux-precice)

add_executable(run_micro_heat_sims EXCLUDE_FROM_ALL micro_sim.cpp)
target_compile_definitions(run_micro_heat_sims PUBLIC "ENABLEMONOLITHIC=0")
target_link_libraries(run_micro_heat_sims PRIVATE dumux-precice)

#add_custom_command(
#    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/micro_sim$(python3 -config --extension-suffix)
#    COMMAND c++ -O3 -Wall -ldumux-precice -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) micro_sim.cpp -o micro_sim$(python3-config --extension-suffix)
#)


#project(micro_sim)

add_subdirectory(pybind11)
pybind11_add_module(micro_sim micro_sim.cpp)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
#target_compile_definitions(micro_sim
#                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
target_link_libraries(micro_sim PRIVATE dumux-precice)

# add a symlink for each input file
add_input_file_links()
