if(CGAL_add_test_included)
  return()
endif(CGAL_add_test_included)
set(CGAL_add_test_included TRUE)

if(NOT POLICY CMP0064)
  # CMake <= 3.3
  if(BUILD_TESTING)
    message(WARNING
      "CGAL CTest support requires CMake 3.4 or later.\n"
      "You must either disable BUILD_TESTING or upgrade CMake.")
  endif()

  # Add a fake function to avoid CMake errors
  function(cgal_add_compilation_test)
  endfunction()

  # Then return, to exit the file
  return()
endif()

include(CTest)

cmake_policy(SET CMP0064 NEW)

include(CMakeParseArguments)

option(CGAL_CTEST_DISPLAY_MEM_AND_TIME
  "Display memory and real time usage at end of CTest test outputs"
  FALSE)
if(CGAL_CTEST_DISPLAY_MEM_AND_TIME)
  find_program(TIME time)
  if(TIME)
    set(TIME_COMMAND ${TIME} -f "\\nMEM: %M\\tTIME: %e\\t%C")
  endif()
endif()

if(ANDROID)
  set(ANDROID_DIR_PREFIX /data/local/tmp/)
  find_program(adb_executable adb)
endif()

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

function(cgal_add_compilation_test exe_name)
  if(NOT POLICY CMP0064 OR TEST compilation_of__${exe_name})
    return()
  endif()
  add_test(NAME "compilation_of__${exe_name}"
    COMMAND ${TIME_COMMAND} "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --target "${exe_name}" --config "$<CONFIG>")
  set_property(TEST "compilation_of__${exe_name}"
    APPEND PROPERTY LABELS "${PROJECT_NAME}")
  if(NOT TARGET ALL_CGAL_TARGETS)
      add_custom_target( ALL_CGAL_TARGETS )
    endif()
  if(NOT TARGET cgal_check_build_system)
    add_custom_target(cgal_check_build_system)
    add_dependencies( ALL_CGAL_TARGETS cgal_check_build_system )
  endif()
  if(NOT TEST check_build_system)
    add_test(NAME "check_build_system"
      COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --target "cgal_check_build_system" --config "$<CONFIG>")
    set_property(TEST "check_build_system"
      APPEND PROPERTY LABELS "Installation")
    if(POLICY CMP0066) # cmake 3.7 or later
      set_property(TEST "check_build_system"
        PROPERTY FIXTURES_SETUP "check_build_system_SetupFixture")
    endif()
  endif()
  if(POLICY CMP0066) # cmake 3.7 or later
    set_property(TEST "compilation_of__${exe_name}"
      APPEND PROPERTY FIXTURES_REQUIRED "check_build_system_SetupFixture")
  endif()
endfunction(cgal_add_compilation_test)

function(cgal_setup_test_properties test_name)
  if(ARGC GREATER 1)
    set(exe_name ${ARGV1})
  endif()
  set_property(TEST "${test_name}"
    APPEND PROPERTY LABELS "${PROJECT_NAME}")
  #      message(STATUS "  working dir: ${CGAL_CURRENT_SOURCE_DIR}")
  set_property(TEST "${test_name}"
    PROPERTY WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  if(exe_name)
    set_property(TEST "${test_name}"
      APPEND PROPERTY DEPENDS "compilation_of__${exe_name}")
  endif()

  get_filename_component(_source_dir_abs ${CMAKE_CURRENT_SOURCE_DIR} ABSOLUTE)
  get_filename_component(_binary_dir_abs ${CMAKE_CURRENT_BINARY_DIR} ABSOLUTE)
  string(FIND "${_binary_dir_abs}" "${_source_dir_abs}" _search_binary_in_source)

  if(_search_binary_in_source EQUAL "-1"
     AND POLICY CMP0066) # CMake 3.7 or later
    if(NOT TEST ${PROJECT_NAME}_SetupFixture)
      if(ANDROID)
        add_test(NAME ${PROJECT_NAME}_SetupFixture
          COMMAND
          ${adb_executable} push
          ${CMAKE_CURRENT_SOURCE_DIR}
          ${ANDROID_DIR_PREFIX}${PROJECT_NAME}
          )
        add_test(NAME ${PROJECT_NAME}_copy_GMP_MPFR
          COMMAND
          ${adb_executable} push
          ${GMP_LIBRARIES} ${MPFR_LIBRARIES}
          ${ANDROID_DIR_PREFIX}${PROJECT_NAME}
          )
        set_property(TEST ${PROJECT_NAME}_copy_GMP_MPFR
          APPEND PROPERTY DEPENDS ${PROJECT_NAME}_SetupFixture)
        set_property(TEST ${PROJECT_NAME}_copy_GMP_MPFR
          PROPERTY FIXTURES_SETUP ${PROJECT_NAME})
      else()
        add_test(NAME ${PROJECT_NAME}_SetupFixture
          COMMAND
          ${CMAKE_COMMAND} -E copy_directory
          ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_BINARY_DIR}/__exec_test_dir
          )
      endif()
      set_property(TEST ${PROJECT_NAME}_SetupFixture
        PROPERTY FIXTURES_SETUP ${PROJECT_NAME})

      if(ANDROID)
        add_test(NAME ${PROJECT_NAME}_CleanupFixture
          COMMAND
          ${adb_executable} shell rm -rf
          ${ANDROID_DIR_PREFIX}${PROJECT_NAME}
          )
      else()
        add_test(NAME ${PROJECT_NAME}_CleanupFixture
          COMMAND
          ${CMAKE_COMMAND} -E remove_directory
          ${CMAKE_CURRENT_BINARY_DIR}/__exec_test_dir
          )
      endif()
      set_property(TEST ${PROJECT_NAME}_CleanupFixture
        PROPERTY FIXTURES_CLEANUP ${PROJECT_NAME})

      set_property(TEST
        ${PROJECT_NAME}_CleanupFixture ${PROJECT_NAME}_SetupFixture
        APPEND PROPERTY LABELS "${PROJECT_NAME}")
    endif()
    if(NOT ANDROID)
      set_property(TEST "${test_name}"
        PROPERTY
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/__exec_test_dir)
    endif()
    set_property(TEST "${test_name}"
      APPEND PROPERTY FIXTURES_REQUIRED "${PROJECT_NAME}")
    if(exe_name)
      set_property(TEST ${test_name}
        APPEND PROPERTY FIXTURES_REQUIRED "${exe_name}")
      set_property(TEST "compilation_of__${exe_name}"
        PROPERTY FIXTURES_SETUP "${exe_name}")
      if(ANDROID)
        add_test(NAME "push_of__${exe_name}"
          COMMAND ${adb_executable} push $<TARGET_FILE:${exe_name}> ${ANDROID_DIR_PREFIX}${PROJECT_NAME}/${exe_name})
        set_property(TEST "push_of__${exe_name}"
          APPEND PROPERTY FIXTURES_SETUP "${exe_name}")
        set_property(TEST "push_of__${exe_name}"
          APPEND PROPERTY DEPENDS "compilation_of__${exe_name}")
      endif()
    endif()
  endif() # end CMake 3.7 or later
endfunction(cgal_setup_test_properties)

function(cgal_add_test exe_name)
  cgal_add_compilation_test(${exe_name})

  cmake_parse_arguments("cgal_add_test" # prefix
    "NO_EXECUTION"                      # optional arguments
    "TEST_NAME"                         # single arguments
    "ARGUMENTS"                         # multivalue arguments
    ${ARGN})
  set(ARGS ${cgal_add_test_ARGUMENTS})
  set(test_name ${cgal_add_test_TEST_NAME})
#  message("  test_name: ${test_name}")
  if(NOT test_name)
    set(test_name "execution___of__${exe_name}")
  endif()
#  message("  test_name: ${test_name}")
  if(cgal_add_test_NO_EXECUTION OR NOT POLICY CMP0064 OR TEST ${test_name})
    return()
  endif()
#  message("Add test ${test_name}")
  set(cin_file "${CGAL_CURRENT_SOURCE_DIR}/${exe_name}.cin")
  if(NOT ARGS AND EXISTS ${cin_file})
    add_test(NAME ${test_name}
      COMMAND ${TIME_COMMAND} ${CMAKE_COMMAND}
      -DCMD:STRING=$<TARGET_FILE:${exe_name}>
      -DCIN:STRING=${cin_file}
      -DANDROID_DIR_PREFIX=${ANDROID_DIR_PREFIX}
      -DPROJECT_NAME=${PROJECT_NAME}
      -P "${CGAL_MODULES_DIR}/run_test_with_cin.cmake")
    set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${cin_file})
    #	message(STATUS "add test: ${exe_name} < ${cin_file}")
  else()
    if(NOT ARGS AND NOT cgal_add_test_TEST_NAME)
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
          string(REGEX REPLACE "\#.*" "" CMD_LINE "${CMD_LINE}")
	  separate_arguments(CMD_LINE_ARGS UNIX_COMMAND ${CMD_LINE})
	  #	    message(STATUS "  args: ${CMD_LINE_ARGS}")
	  list(APPEND ARGS ${CMD_LINE_ARGS})
        endforeach()
        expand_list_with_globbing(ARGS)
        set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${cmd_file})
      endif()
    endif()
    #	message(STATUS "add test: ${exe_name} ${ARGS}")
    if(ANDROID)
      add_test(NAME ${test_name} COMMAND ${TIME_COMMAND} ${adb_executable} shell cd ${ANDROID_DIR_PREFIX}${PROJECT_NAME} && LD_LIBRARY_PATH=${ANDROID_DIR_PREFIX}${PROJECT_NAME} ${ANDROID_DIR_PREFIX}${PROJECT_NAME}/${exe_name} ${ARGS})
    else()
      add_test(NAME ${test_name} COMMAND ${TIME_COMMAND} $<TARGET_FILE:${exe_name}> ${ARGS})
    endif()
  endif()
  cgal_setup_test_properties(${test_name} ${exe_name})
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
  add_test(NAME ${test_name}
    COMMAND ${TIME_COMMAND} $<TARGET_FILE:${exe_name}> ${ARGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set_property(TEST "${test_name}"
    APPEND PROPERTY LABELS "${PROJECT_NAME}")
endfunction()

function(CGAL_add_compilation_tests_for_all_targets)
endfunction()
