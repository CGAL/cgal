if(CGAL_add_test_included)
  return()
endif(CGAL_add_test_included)
set(CGAL_add_test_included TRUE)

# Enable testing with CGAL_ENABLE_TESTING. Before CGAL-6.0, users would enable
# the tests by specifying BUILD_TESTING. For compatibility, If BUILD_TESTING is
# set, that is the default value for CGAL_ENABLE_TESTING. Otherwise, the default
# value is OFF.
option(CGAL_ENABLE_TESTING "Build the testing tree." ${BUILD_TESTING})

if(CGAL_ENABLE_TESTING)
  enable_testing()
endif()


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
  set(CGAL_REMOTE_TEST_DIR_PREFIX /data/local/tmp/ CACHE PATH "Path to the directory where the tests will be executed in a remote testsuite.")
  find_program(adb_executable adb)
endif()
if(CGAL_RUN_TESTS_THROUGH_SSH)
  set(CGAL_REMOTE_TEST_DIR_PREFIX /home/pi/CGAL/ CACHE PATH "Path to the directory where the tests will be executed in a remote testsuite.")
  find_program(ssh_executable ssh)
  find_program(scp_executable scp)
endif()

# Process a list, and replace items contains a file pattern (like
# `*.off`) by the sublist that corresponds to the globbing of the
# pattern in the directory `${CMAKE_CURRENT_SOURCE_DIR}`.
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
      file(GLOB files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${item_n})
      list(APPEND output_list ${files})
    else()
      list(APPEND output_list ${item_n})
    endif()
#    message(STATUS "  new value of the output list: ${output_list}")
  endforeach()
  set(${list_name} ${output_list} PARENT_SCOPE)
endfunction()

function(cgal_add_compilation_test exe_name)
  cmake_policy(SET CMP0064 NEW)
  if(NOT CMAKE_VS_MSBUILD_COMMAND)
    if(TEST "compilation of  ${exe_name}")
      return()
    endif()
    add_test(NAME "compilation of  ${exe_name}"
      COMMAND ${TIME_COMMAND} "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --target "${exe_name}" --config "$<CONFIG>")
    set_property(TEST "compilation of  ${exe_name}"
      APPEND PROPERTY LABELS "${PROJECT_NAME}")
    set_property(TEST "compilation of  ${exe_name}"
      APPEND PROPERTY FIXTURES_REQUIRED "check_build_system_SetupFixture")
  elseif(NOT TARGET "compilation_of__${PROJECT_NAME}")#CMAKE_VS_MSBUILD_COMMAND
    #this target is just a flag, to deal with the scope problem with the tests
    add_custom_target("compilation_of__${PROJECT_NAME}")
    add_test(NAME "compilation of  ${PROJECT_NAME}"
      COMMAND ${TIME_COMMAND} "${CMAKE_VS_MSBUILD_COMMAND}" "${PROJECT_BINARY_DIR}/${PROJECT_NAME}.sln" "-m:$ENV{NUMBER_OF_PROCESSORS}" "/t:Build" "/p:Configuration=$<CONFIG>")
    set_property(TEST "compilation of  ${PROJECT_NAME}"
      APPEND PROPERTY LABELS "${PROJECT_NAME}")
    set_property(TEST "compilation of  ${PROJECT_NAME}"
      APPEND PROPERTY FIXTURES_REQUIRED "check_build_system_SetupFixture")
    set_tests_properties("compilation of  ${PROJECT_NAME}"
      PROPERTIES RUN_SERIAL TRUE)
    #because of the scope of the tests, this part cannot go in the relevant CMakeLists
    if("${PROJECT_NAME}" STREQUAL "Polyhedron_Demo")
      set_tests_properties("compilation of  Polyhedron_Demo" PROPERTIES TIMEOUT 2400)
    elseif("${PROJECT_NAME}" STREQUAL "Mesh_3_Tests" OR "${PROJECT_NAME}" STREQUAL "Mesh_3_Examples")
      set_tests_properties("compilation of  ${PROJECT_NAME}" PROPERTIES TIMEOUT 1600)
    endif()

  endif()#CMAKE_VS_MSBUILD_COMMAND
  if(NOT TARGET ALL_CGAL_TARGETS)
    add_custom_target( ALL_CGAL_TARGETS )
  endif()
  if(NOT TARGET cgal_check_build_system)
    add_custom_target(cgal_check_build_system)
    add_dependencies( ALL_CGAL_TARGETS cgal_check_build_system )
    add_test(NAME "check build system"
      COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --target "cgal_check_build_system" --config "$<CONFIG>")
    set_property(TEST "check build system"
      APPEND PROPERTY LABELS "CGAL_build_system" "Installation")
    set_property(TEST "check build system"
       PROPERTY FIXTURES_SETUP "check_build_system_SetupFixture")
  endif()
  if(TARGET CGAL_Qt6_moc_and_resources) # if CGAL_Qt6 was searched, and is header-only
    get_property(linked_libraries TARGET "${exe_name}" PROPERTY LINK_LIBRARIES)
    #  message(STATUS "${exe_name} depends on ${linked_libraries}")
    if( ("CGAL_Qt6" IN_LIST linked_libraries OR "CGAL::CGAL_Qt6" IN_LIST linked_libraries OR "CGAL::CGAL_Basic_viewer" IN_LIST linked_libraries)
        AND NOT TARGET "compilation_of__CGAL_Qt6_moc_and_resources")
      # This custom target is useless. It is used only as a flag to
      # detect that the test has already been created.
      add_custom_target("compilation_of__CGAL_Qt6_moc_and_resources")
      add_dependencies( "compilation_of__CGAL_Qt6_moc_and_resources" CGAL_Qt6_moc_and_resources )
      add_test(NAME "compilation of  CGAL_Qt6_moc_and_resources"
        COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --target "CGAL_Qt6_moc_and_resources" --config "$<CONFIG>")
      set_property(TEST "compilation of  CGAL_Qt6_moc_and_resources"
        APPEND PROPERTY LABELS "CGAL_build_system" "Installation")
      set_property(TEST "compilation of  CGAL_Qt6_moc_and_resources"
        PROPERTY FIXTURES_SETUP "CGAL_Qt6_moc_and_resources_Fixture")
      set_property(TEST "compilation of  CGAL_Qt6_moc_and_resources"
        APPEND PROPERTY DEPENDS "check build system")
    endif()
  endif()
endfunction(cgal_add_compilation_test)

option(CGAL_TEST_DRAW_FUNCTIONS "If set, the ctest command will not skip the tests of the draw functions.")

function(cgal_setup_test_properties test_name)
  if(ARGC GREATER 1)
    set(exe_name ${ARGV1})
  endif()
  set_property(TEST "${test_name}"
    APPEND PROPERTY LABELS "${PROJECT_NAME}")
  #      message(STATUS "  working dir: ${CMAKE_CURRENT_SOURCE_DIR}")
  set_property(TEST "${test_name}"
    PROPERTY WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  if(NOT CGAL_TEST_DRAW_FUNCTIONS)
    set_property(TEST "${test_name}"
      APPEND PROPERTY ENVIRONMENT CGAL_TEST_SUITE=1)
  endif()

  if(exe_name)
    if(NOT CMAKE_VS_MSBUILD_COMMAND)
      set_property(TEST "${test_name}"
        APPEND PROPERTY DEPENDS "compilation of  ${exe_name}")
    else()#CMAKE_VS_MSBUILD_COMMAND
      set_property(TEST "${test_name}"
        APPEND PROPERTY DEPENDS "compilation of  ${PROJECT_NAME}")
    endif()#CMAKE_VS_MSBUILD_COMMAND
  endif()

  get_filename_component(_source_dir_abs ${CMAKE_CURRENT_SOURCE_DIR} ABSOLUTE)
  get_filename_component(_binary_dir_abs ${CMAKE_CURRENT_BINARY_DIR} ABSOLUTE)
  string(FIND "${_binary_dir_abs}" "${_source_dir_abs}" _search_binary_in_source)

  if(_search_binary_in_source EQUAL "-1")
    if(NOT TEST "copy source_dir of ${PROJECT_NAME}")
      if(ANDROID)
        add_test(NAME "copy source_dir of ${PROJECT_NAME}"
          COMMAND
          ${adb_executable} push
          ${CMAKE_CURRENT_SOURCE_DIR}
          ${CGAL_REMOTE_TEST_DIR_PREFIX}${PROJECT_NAME}
          )
        add_test(NAME ${PROJECT_NAME}_copy_GMP_MPFR
          COMMAND
          ${adb_executable} push
          ${GMP_LIBRARIES} ${MPFR_LIBRARIES}
          ${CGAL_REMOTE_TEST_DIR_PREFIX}${PROJECT_NAME}
          )
        set_property(TEST ${PROJECT_NAME}_copy_GMP_MPFR
          APPEND PROPERTY DEPENDS "copy source_dir of ${PROJECT_NAME}")
        set_property(TEST ${PROJECT_NAME}_copy_GMP_MPFR
          PROPERTY FIXTURES_SETUP ${PROJECT_NAME})
      elseif(CGAL_RUN_TESTS_THROUGH_SSH)
        add_test(NAME "copy source_dir of ${PROJECT_NAME}"
          COMMAND
          ${scp_executable} -r
          ${CMAKE_CURRENT_SOURCE_DIR}
          ${SSH_HOST}:${CGAL_REMOTE_TEST_DIR_PREFIX}${PROJECT_NAME}
          )
      else()
        add_test(NAME "copy source_dir of ${PROJECT_NAME}"
          COMMAND
          ${CMAKE_COMMAND} -E copy_directory
          ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_BINARY_DIR}/__exec_test_dir
          )
      endif()
      set_property(TEST "copy source_dir of ${PROJECT_NAME}"
        PROPERTY FIXTURES_SETUP ${PROJECT_NAME})

      if(ANDROID)
        add_test(NAME "cleanup of ${PROJECT_NAME}"
          COMMAND
          ${adb_executable} shell rm -rf
          ${CGAL_REMOTE_TEST_DIR_PREFIX}${PROJECT_NAME}
          )
      elseif(CGAL_RUN_TESTS_THROUGH_SSH)
        add_test(NAME "cleanup of ${PROJECT_NAME}"
          COMMAND
          ${ssh_executable} ${SSH_HOST} rm -rf
          ${CGAL_REMOTE_TEST_DIR_PREFIX}${PROJECT_NAME}
          )
      else()
        add_test(NAME "cleanup of ${PROJECT_NAME}"
          COMMAND
          ${CMAKE_COMMAND} -E remove_directory
          ${CMAKE_CURRENT_BINARY_DIR}/__exec_test_dir
          )
      endif()
      set_property(TEST "cleanup of ${PROJECT_NAME}"
        PROPERTY FIXTURES_CLEANUP ${PROJECT_NAME})

      set_property(TEST
        "cleanup of ${PROJECT_NAME}" "copy source_dir of ${PROJECT_NAME}"
        APPEND PROPERTY LABELS "${PROJECT_NAME}")
    endif()
    if(NOT ANDROID AND NOT CGAL_RUN_TESTS_THROUGH_SSH)
      set_property(TEST "${test_name}"
        PROPERTY
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/__exec_test_dir)
    endif()

    set_property(TEST "${test_name}"
      APPEND PROPERTY FIXTURES_REQUIRED "${PROJECT_NAME}")
    if(exe_name)
      set_property(TEST ${test_name}
        APPEND PROPERTY FIXTURES_REQUIRED "${exe_name}")
      if(NOT CMAKE_VS_MSBUILD_COMMAND)
        set_property(TEST "compilation of  ${exe_name}"
          PROPERTY FIXTURES_SETUP "${exe_name}")
      else()#CMAKE_VS_MSBUILD_COMMAND
        set_property(TEST "compilation of  ${PROJECT_NAME}"
          PROPERTY FIXTURES_SETUP "${exe_name}")
      endif()#CMAKE_VS_MSBUILD_COMMAND
      if((ANDROID OR CGAL_RUN_TESTS_THROUGH_SSH) AND NOT TEST push_of__${exe_name})
        if(ANDROID)
          add_test(NAME "push_of__${exe_name}"
            COMMAND ${adb_executable} push $<TARGET_FILE:${exe_name}> ${CGAL_REMOTE_TEST_DIR_PREFIX}${PROJECT_NAME}/${exe_name})
        elseif(CGAL_RUN_TESTS_THROUGH_SSH)
          add_test(NAME "push_of__${exe_name}"
            COMMAND ${scp_executable} $<TARGET_FILE:${exe_name}> ${SSH_HOST}:${CGAL_REMOTE_TEST_DIR_PREFIX}${PROJECT_NAME}/)
        endif()
        set_property(TEST "push_of__${exe_name}"
          APPEND PROPERTY DEPENDS "compilation of  ${exe_name}")
        set_property(TEST "push_of__${exe_name}"
          APPEND PROPERTY FIXTURES_SETUP "${exe_name}")
        set_property(TEST "push_of__${exe_name}"
          APPEND PROPERTY FIXTURES_REQUIRED "${PROJECT_NAME}")
        set_property(TEST "push_of__${exe_name}"
          APPEND PROPERTY LABELS "${PROJECT_NAME}")
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
    set(test_name "execution   of  ${exe_name}")
  endif()
#  message("  test_name: ${test_name}")
  if(cgal_add_test_NO_EXECUTION OR TEST ${test_name})
    return()
  endif()
#  message("Add test ${test_name}")
  set(cin_file "${CMAKE_CURRENT_SOURCE_DIR}/${exe_name}.cin")
  if(NOT ARGS AND EXISTS ${cin_file})
    if(CGAL_RUN_TESTS_THROUGH_SSH)
      add_test(NAME ${test_name} COMMAND bash -c "${TIME_COMMAND} ${ssh_executable}  ${SSH_HOST} \"cd ${CGAL_REMOTE_TEST_DIR_PREFIX}${PROJECT_NAME} && ${CGAL_REMOTE_TEST_DIR_PREFIX}${PROJECT_NAME}/${exe_name} 3< <(cat; kill -INT 0) < ${exe_name}.cin\" <&1")
    else()
      if(ANDROID)
        set(cmd ${exe_name})
      else()
        set(cmd $<TARGET_FILE:${exe_name}>)
      endif()
      add_test(NAME ${test_name}
        COMMAND ${TIME_COMMAND} ${CMAKE_COMMAND}
        -DCMD:STRING=${cmd}
        -DCIN:STRING=${cin_file}
        -DCGAL_REMOTE_TEST_DIR_PREFIX=${CGAL_REMOTE_TEST_DIR_PREFIX}
        -DCGAL_RUN_TESTS_THROUGH_SSH=${CGAL_RUN_TESTS_THROUGH_SSH}
        -DSSH_HOST=${SSH_HOST}
        -Dssh_executable=${ssh_executable}
        -DCGAL_REMOTE_TEST_DIR_PREFIX=${CGAL_REMOTE_TEST_DIR_PREFIX}
        -DPROJECT_NAME=${PROJECT_NAME}
        -P "${CGAL_MODULES_DIR}/run_test_with_cin.cmake")
      set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${cin_file})
      # message(STATUS "add test: ${exe_name} < ${cin_file}")
    endif()
  else()
    if(NOT ARGS AND NOT cgal_add_test_TEST_NAME)
      if(ARGC GREATER 2 AND ARGV2)
        set(cmd_file "${CMAKE_CURRENT_SOURCE_DIR}/${ARGV2}.cmd")
      elseif(ARGC GREATER 1 AND ARGV1 AND NOT EXISTS ${cmd_file})
        set(cmd_file "${CMAKE_CURRENT_SOURCE_DIR}/${ARGV1}.cmd")
      elseif(NOT EXISTS ${cmd_file})
        set(cmd_file "${CMAKE_CURRENT_SOURCE_DIR}/${exe_name}.cmd")
      endif()
      if(EXISTS ${cmd_file})
        file(STRINGS "${cmd_file}" CMD_LINES)
        string(CONFIGURE "${CMD_LINES}" CMD_LINES)
        set(ARGS)
        #         message(STATUS "DEBUG test ${exe_name}")
        foreach(CMD_LINE ${CMD_LINES})
    #       message(STATUS "  command line: ${CMD_LINE}")
          string(REGEX REPLACE "\#.*" "" CMD_LINE "${CMD_LINE}")
          separate_arguments(CMD_LINE_ARGS UNIX_COMMAND ${CMD_LINE})
    #       message(STATUS "  args: ${CMD_LINE_ARGS}")
          list(APPEND ARGS ${CMD_LINE_ARGS})
        endforeach()
        expand_list_with_globbing(ARGS)
        set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${cmd_file})
      endif()
    endif()
    #   message(STATUS "add test: ${exe_name} ${ARGS}")
    if(ANDROID)
      add_test(NAME ${test_name} COMMAND ${TIME_COMMAND} ${adb_executable} shell cd ${CGAL_REMOTE_TEST_DIR_PREFIX}${PROJECT_NAME} && LD_LIBRARY_PATH=${CGAL_REMOTE_TEST_DIR_PREFIX}${PROJECT_NAME} ${CGAL_REMOTE_TEST_DIR_PREFIX}${PROJECT_NAME}/${exe_name} ${ARGS})
    elseif(CGAL_RUN_TESTS_THROUGH_SSH)
      STRING(REPLACE ";" " " arg_str "${ARGS}")
      add_test(NAME ${test_name} COMMAND bash -c "${TIME_COMMAND} ${ssh_executable}  ${SSH_HOST} \"cd ${CGAL_REMOTE_TEST_DIR_PREFIX}${PROJECT_NAME} && ${CGAL_REMOTE_TEST_DIR_PREFIX}${PROJECT_NAME}/${exe_name} ${arg_str} 3< <(cat; kill -INT 0)\" <&1")
    else()
      add_test(NAME ${test_name} COMMAND ${TIME_COMMAND} $<TARGET_FILE:${exe_name}> ${ARGS})
    endif()
  endif()
  cgal_setup_test_properties(${test_name} ${exe_name})
endfunction()

function(CGAL_add_compilation_tests_for_all_targets)
endfunction()
