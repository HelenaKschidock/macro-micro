add_executable(test_micro_heat EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_micro_heat PUBLIC "ENABLEMONOLITHIC=0")
target_link_libraries(test_micro_heat PRIVATE dumux-precice)

# add a symlink for each input file
add_input_file_links()
