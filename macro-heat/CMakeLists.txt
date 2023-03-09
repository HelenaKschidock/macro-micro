add_executable(test_macro_heat EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_macro_heat PUBLIC "ENABLEMONOLITHIC=0")
target_link_libraries(test_macro_heat PRIVATE dumux-precice)
target_compile_definitions(test_macro_heat PUBLIC TYPETAG=OnePNIConductionCCTpfa)

# add a symlink for each input file
add_input_file_links()

macro(add_precice_file_links)
  FILE(GLOB precice_input_files *.xml)
  foreach(VAR ${precice_input_files})
    get_filename_component(file_name ${VAR} NAME)
    dune_symlink_to_source_files(FILES ${file_name})
  endforeach()
endmacro()
add_precice_file_links()
