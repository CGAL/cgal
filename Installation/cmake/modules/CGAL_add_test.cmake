if(CGAL_add_test_included)
  return()
endif(CGAL_add_test_included)
set(CGAL_add_test_included TRUE)


# Process a list, and replace items contains a file pattern (like
# `*.off`) by the sublist that corresponds to the globbing of the
# pattern in the directory `${CGAL_CURRENT_SOURCE_DIR}`.
#
#
# For example: the `file
# test/Poisson_surface_reconstruction_3/poisson_reconstruction_test.cmd`
# contains:
#
#     data/*.off data/*.xyz data/*.pwn
#
# For that file, the list `ARGS` computed in the function
# `create_single_source_cgal_program` (see below) is the list:
#
#     data/*.off;data/*.xyz;data/*.pwn
#
# A call to `expand_list_with_globbing(ARGS)` replaces the list by:
#
#     data/ChineseDragon-10kv.off;data/robocat_deci.off;data/sphere_20k.xyz;data/oni.pwn"
#
function(expand_list_with_globbing list_name)
  set(input_list ${${list_name}})
#  message(STATUS "expand_list_with_globbing(${list_name}), ${list_name} is: ${input_list}")
  list(LENGTH input_list list_lenght)
  math(EXPR list_last_n "${list_lenght} - 1")
  set(output_list)
  foreach(n RANGE ${list_last_n})
#    message(STATUS "n=${n}")
    list(GET input_list ${n} item_n)
#    message(STATUS "argument ${n} is ${item_n}")
    if(item_n MATCHES ".*\\*.*")
      file(GLOB files RELATIVE ${CGAL_CURRENT_SOURCE_DIR} ${item_n})
      list(APPEND output_list ${files})
    else()
      list(APPEND output_list ${item_n})
    endif()
#    message(STATUS "  new value of the output list: ${output_list}")
  endforeach()
  set(${list_name} ${output_list} PARENT_SCOPE)
endfunction()

function(cgal_add_test exe_name)
  set(cin_file "${CGAL_CURRENT_SOURCE_DIR}/${exe_name}.cin")
  if(EXISTS ${cin_file})
    add_test(NAME ${exe_name}
      COMMAND ${CMAKE_COMMAND}
      -DCMD:STRING=$<TARGET_FILE:${exe_name}>
      -DCIN:STRING=${cin_file}
      -P "${CGAL_MODULES_DIR}/run_test_with_cin.cmake")
    #	message(STATUS "add test: ${exe_name} < ${cin_file}")
  else()
    # TODO: deal with shell globbing; if the `cmd` file contains
    # a `*`, then interprete the command using bash.
    if(ARGC GREATER 2 AND ARGV2)
      set(cmd_file "${CGAL_CURRENT_SOURCE_DIR}/${ARGV2}.cmd")
    elseif(ARGC GREATER 1 AND ARGV1 AND NOT EXISTS ${cmd_file})
      set(cmd_file "${CGAL_CURRENT_SOURCE_DIR}/${ARGV1}.cmd")
    elseif(NOT EXISTS ${cmd_file})
      set(cmd_file "${CGAL_CURRENT_SOURCE_DIR}/${exe_name}.cmd")
    endif()
    if(EXISTS ${cmd_file})
      file(STRINGS "${cmd_file}" CMD_LINES)
      set(ARGS)
      #	  message(STATUS "DEBUG test ${exe_name}")
      foreach(CMD_LINE ${CMD_LINES})
	#	    message(STATUS "  command line: ${CMD_LINE}")
	separate_arguments(CMD_LINE_ARGS UNIX_COMMAND ${CMD_LINE})
	#	    message(STATUS "  args: ${CMD_LINE_ARGS}")
	list(APPEND ARGS ${CMD_LINE_ARGS})
      endforeach()
      expand_list_with_globbing(ARGS)
    endif()
    #	message(STATUS "add test: ${exe_name} ${ARGS}")
    add_test(NAME ${exe_name} COMMAND ${exe_name} ${ARGS})
  endif()
  set_property(TEST "${exe_name}"
    APPEND PROPERTY LABELS "${PROJECT_NAME}")
  #      message(STATUS "  working dir: ${CGAL_CURRENT_SOURCE_DIR}")
  set_property(TEST "${exe_name}"
    PROPERTY WORKING_DIRECTORY ${CGAL_CURRENT_SOURCE_DIR})
  add_test(NAME "build_target_${exe_name}"
    COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --target "${exe_name}")
  set_property(TEST "build_target_${exe_name}"
    APPEND PROPERTY LABELS "${PROJECT_NAME}")
  set_property(TEST "${exe_name}"
    APPEND PROPERTY DEPENDS "build_target_${exe_name}")

  return()

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${exe_name}.cin")
    set(ARGS "${CMAKE_CURRENT_SOURCE_DIR}/${exe_name}.cin")
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${ARGV2}.cmd")
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${ARGV2}.cmd"
      ARGS LIMIT_COUNT 1)
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${ARGV1}.cmd")
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${ARGV1}.cmd"
      ARGS LIMIT_COUNT 1)
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${exe_name}.cmd")
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/${exe_name}.cmd"
      ARGS LIMIT_COUNT 1)
    # TODO: handle multi-lines .cmd files
    # see https://github.com/CGAL/cgal/pull/1295/files/c65d3abe17bb3e677b8077996cdaf8672f9c4c6f#r71705451
  endif()
  string(REPLACE ";" " " args_str "${ARGS}")
  add_test(NAME ${exe_name}
    COMMAND ${exe_name} ${ARGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set_property(TEST "${exe_name}"
    APPEND PROPERTY LABELS "${PROJECT_NAME}")
endfunction()
