add_executable(test_micro_heat EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_micro_heat PUBLIC "ENABLEMONOLITHIC=0")
target_link_libraries(test_micro_heat PRIVATE dumux-precice)

add_executable(run_micro_heat_sims EXCLUDE_FROM_ALL micro_sim.cpp)
target_compile_definitions(run_micro_heat_sims PUBLIC "ENABLEMONOLITHIC=0")
target_link_libraries(run_micro_heat_sims PRIVATE dumux-precice)

add_custom_target(copy_micro_heat_all ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

#macro(add_micro_manager_file_links)
#    # create a list of files to copy
#    set( micro_manager_files
#        ${CMAKE_CURRENT_SOURCE_DIR}/micro-manager-config.json
#        ${CMAKE_CURRENT_SOURCE_DIR}/run_micro_manager.py 
#    )
#  foreach(VAR ${micro_manager_files})
#    get_filename_component(file_name ${VAR} NAME)
#    dune_symlink_to_source_files(FILES ${file_name})
#  endforeach()
#endmacro()
#add_micro_manager_file_links()

# add a symlink for each input file
add_input_file_links()
