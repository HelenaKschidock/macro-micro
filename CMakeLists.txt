add_subdirectory(macro-heat)
add_subdirectory(nutils-micro-heat)
add_subdirectory(micro-heat)

add_custom_target(copy_jobscript ALL
    COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/macro_micro_jobscript.sh
            ${CMAKE_CURRENT_BINARY_DIR}/macro_micro_jobscript.sh)


