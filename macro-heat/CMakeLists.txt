add_executable(test_macro_heat EXCLUDE_FROM_ALL main_macro_heat.cc)
target_compile_definitions(test_macro_heat PUBLIC "ENABLEMONOLITHIC=0")
target_link_libraries(test_macro_heat PRIVATE dumux-precice)

# add a symlink for each input file
add_input_file_links()

macro(add_precice_file_links)
  FILE(GLOB precice_files *.xml)
  foreach(VAR ${precice_input_files})
    get_filename_component(file_name ${VAR} NAME)
    dune_symlink_to_source_files(FILES ${file_name})
  endforeach()
endmacro()
add_precice_file_links()

# Dummy test is needed to build porous media solver. The default
# `dumux_add_test` (and the underlying `dune_add_test`) only accept
# one target executable.
dumux_add_test(NAME dummy_test_macro_heat
              TARGET test_macro_heat
              TIMEOUT 30
              COMMAND ${CMAKE_SOURCE_DIR}/test/return-test-passed.sh
)

dumux_add_test(NAME test_macro_heat
              TARGET test_macro_heat
              TIMEOUT 60
              COMMAND ${CMAKE_SOURCE_DIR}/test/run-iterative-test.py
              CMD_ARGS
              --dumux-param-file ${CMAKE_SOURCE_DIR}/examples/macro-micro/params.input 
              --precice-config-file ${CMAKE_SOURCE_DIR}/examples/macro-micro/precice-config.xml
)  