cmake_minimum_required(VERSION 3.19...3.31)

project( Installation_Tests )


macro(create_link_to_program COMPONENT)
  add_executable(link_to_${COMPONENT} link_to_${COMPONENT}.cpp)
  target_link_libraries(link_to_${COMPONENT} CGAL::${COMPONENT})
  add_to_cached_list(CGAL_EXECUTABLE_TARGETS link_to_${COMPONENT})
  CGAL_add_test(link_to_${COMPONENT})
endmacro()

find_package(CGAL REQUIRED COMPONENTS Core)

include(CGAL_CreateSingleSourceCGALProgram)
include(CGAL_SetupGMP)

create_single_source_cgal_program("endian.cpp")

create_single_source_cgal_program("deprecation_warning.cpp")

create_single_source_cgal_program("test_use_h.cpp")

create_single_source_cgal_program(will_fail.cpp will_fail_aux.cpp)
if(CGAL_ENABLE_TESTING)
  set_property(TEST "execution   of  will_fail" PROPERTY WILL_FAIL TRUE)
endif()

find_package(TBB QUIET)
include(CGAL_TBB_support)
create_single_source_cgal_program("test_TBB.cpp")
if(TARGET CGAL::TBB_support)
  target_link_libraries(test_TBB PRIVATE CGAL::TBB_support)
endif()

create_link_to_program(CGAL)

if ( CGAL_Core_FOUND )
  create_link_to_program(CGAL_Core)
endif()

if(WITH_CGAL_ImageIO)
  create_link_to_program(CGAL_ImageIO)
endif()


if(WIN32 OR CMAKE_SYSTEM_NAME STREQUAL Windows)
  create_single_source_cgal_program(test_gmp_mpfr_dll.cpp)
  target_link_libraries(test_gmp_mpfr_dll PRIVATE version)
  string(RANDOM RDM_DIR)#5 random chars to avoid conflicts in parallel testsuites
  set(WORKING_DIRS_PREFIX $ENV{TEMP}/${RDM_DIR})
else()
  set(WORKING_DIRS_PREFIX ${CMAKE_BINARY_DIR})
endif()

if(WITH_CGAL_Qt6)
  find_package(Qt6 QUIET COMPONENTS Core)
  if(Qt6_FOUND)
    create_link_to_program(CGAL_Qt6)
  endif()
endif()

find_package(LEDA QUIET)
if(LEDA_FOUND)
  include(${LEDA_USE_FILE})
  add_executable(test_LEDA_IO test_LEDA_IO.cpp)
  target_link_libraries(test_LEDA_IO ${LEDA_LIBRARIES})
endif()

if("$ENV{INIT_FILE}" STREQUAL "")
  set(INIT_FILE)
else()
  set(INIT_FILE "-C$ENV{INIT_FILE}")
endif()
##
# Function: CGAL_installation_test_find_package_version
# Adds a test for finding the CGAL package with a specific version mode.
##
function(CGAL_installation_test_find_package_version mode)
  file(REMOVE_RECURSE "${CMAKE_CURRENT_BINARY_DIR}/build-test_find_package_version_${mode}")
  file(REMOVE_RECURSE "${CMAKE_CURRENT_BINARY_DIR}/test_find_package_version_${mode}")
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/build-test_find_package_version_${mode}")
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_find_package_version_${mode}")

  set(EXACT "")
  set(FAIL FALSE)
  if(mode STREQUAL "less")
    set(CGAL_MINOR_VERSION 0)
  elseif(mode STREQUAL "less_major")
    set(CGAL_MAJOR_VERSION 0)
  elseif(mode STREQUAL "greater" OR mode STREQUAL "fail-exact")
    math(EXPR CGAL_MINOR_VERSION "${CGAL_MINOR_VERSION} + 1")
  endif()
  if(mode STREQUAL "exact" OR mode STREQUAL "fail-exact")
    set(EXACT "EXACT ")
  endif()
  if(ARGN STREQUAL "FAIL")
    set(FAIL TRUE)
  endif()

  set(VERSION "${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}")
  if(CGAL_BUGFIX_VERSION AND CGAL_BUGFIX_VERSION GREATER 0)
    set(VERSION "${VERSION}.${CGAL_BUGFIX_VERSION}")
  endif()
  configure_file(test_find_package.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/test_find_package_version_${mode}/CMakeLists.txt" @ONLY)

  add_test(NAME "test_find_package_version_${mode}"
    COMMAND ${CMAKE_COMMAND} ${INIT_FILE}
                             -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                             -S "${CMAKE_CURRENT_BINARY_DIR}/test_find_package_version_${mode}"
                             -B "${CMAKE_CURRENT_BINARY_DIR}/build-test_find_package_version_${mode}"
    )
  cgal_setup_test_properties("test_find_package_version_${mode}")
  if(ARGN STREQUAL "FAIL")
    set_tests_properties("test_find_package_version_${mode}" PROPERTIES WILL_FAIL TRUE)
  endif()
endfunction()

CGAL_installation_test_find_package_version(less)
CGAL_installation_test_find_package_version(less_major)
CGAL_installation_test_find_package_version(equal)
CGAL_installation_test_find_package_version(greater FAIL)
CGAL_installation_test_find_package_version(exact)
CGAL_installation_test_find_package_version(fail-exact FAIL)

find_package(CGAL COMPONENTS Qt6)

#If ctest is ran from a global config, CGAL_SOURCE_DIR exists, but from Installation/test it doesn't. In that case, however, there is a CGAL_DIR.
if(NOT CGAL_SOURCE_DIR)
  if(NOT ENV{CGAL_RELEASE_DIR})
    set(CGAL_SOURCE_DIR "${CGAL_DIR}")
  else()
    set(CGAL_SOURCE_DIR "${ENV{CGAL_RELEASE_DIR}}")
  endif()
endif()
string(REGEX REPLACE "/cygdrive/c/" "C:/" CGAL_SOURCE_DIR "${CGAL_SOURCE_DIR}")

file(MAKE_DIRECTORY "${WORKING_DIRS_PREFIX}/test_config_file")
configure_file(test_configuration.cmake.in "${CMAKE_BINARY_DIR}/test_config_file/CMakeLists.txt" @ONLY)

if(CGAL_Qt6_FOUND)
  file(MAKE_DIRECTORY "${WORKING_DIRS_PREFIX}/test_config_file_qt")
  configure_file(test_configuration_qt.cmake.in "${CMAKE_BINARY_DIR}/test_config_file_qt/CMakeLists.txt" @ONLY)
endif()#CGAL_Qt6_FOUND

if(EXISTS "${CGAL_SOURCE_DIR}/Installation/lib/cmake/CGAL/CGALConfig.cmake")
  set(CGAL_INSTALLATION_FILES_DIR "${CGAL_SOURCE_DIR}/Installation")
else()
  set(CGAL_INSTALLATION_FILES_DIR "${CGAL_SOURCE_DIR}")
endif()

##
# Function: CGAL_add_config_file_test
# Adds a config file test for a given CGAL config path
##
function(CGAL_add_config_file_test test_name src_dir cgal_config_path)
  cmake_parse_arguments(PARSE_ARGV 3 arg "TEST_EXAMPLES" "" "FIXTURES_REQUIRED" )
  if(arg_FIXTURES_REQUIRED)
    set(_recursive_call_extra_args FIXTURES_REQUIRED ${arg_FIXTURES_REQUIRED})
  else()
    set(_recursive_call_extra_args)
  endif()
  list(APPEND arg_FIXTURES_REQUIRED "${test_name}_fixture")
  set(build_dir "${CMAKE_BINARY_DIR}/CGAL_cmake_test_suite/build-${test_name}")
  get_filename_component(cgal_dir "${cgal_config_path}" DIRECTORY)
  if(arg_TEST_EXAMPLES)
    set(extra_ctest_args -R terrain_with_info)
  endif()
  add_test(NAME "config___${test_name}"
    COMMAND ${CMAKE_COMMAND} ${GENERATOR} ${INIT_FILE}
                             -DCGAL_ENABLE_TESTING=ON
                             -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                             "-DCGAL_DIR=${cgal_dir}"
                             "-DCGAL_GIVEN_DIR=${cgal_dir}"
                             -S "${src_dir}" -B "${build_dir}")
  add_test(NAME "setup____${test_name}"  COMMAND ${CMAKE_COMMAND} -E make_directory "${build_dir}")
  add_test(NAME "cleanup__${test_name}"  COMMAND ${CMAKE_COMMAND} -E remove_directory "${build_dir}")
  add_test(NAME "test_____${test_name}"
    COMMAND ${CMAKE_CTEST_COMMAND} --test-dir "${build_dir}"
                                   --build-config "$<CONFIG>"
                                   --output-on-failure
                                   ${extra_ctest_args})
  set_property(TEST "test_____${test_name}" APPEND PROPERTY DEPENDS "config___${test_name}")
  set_property(TEST "test_____${test_name}" APPEND PROPERTY ENVIRONMENT CTEST_NO_TESTS_ACTION=error)

  set_property(TEST "setup____${test_name}" APPEND PROPERTY FIXTURES_SETUP "${test_name}_fixture")
  set_property(TEST "cleanup__${test_name}" APPEND PROPERTY FIXTURES_CLEANUP "${test_name}_fixture")
  set_property(TEST "config___${test_name}"
                    "test_____${test_name}"
              APPEND PROPERTY FIXTURES_REQUIRED ${arg_FIXTURES_REQUIRED})

  set_property(TEST "setup____${test_name}"
                    "cleanup__${test_name}"
                    "config___${test_name}"
                    "test_____${test_name}"
               APPEND PROPERTY LABELS Installation_Tests CGAL_cmake_testsuite)

  if(NOT arg_TEST_EXAMPLES)
    # recursive call to test examples
    CGAL_add_config_file_test(
      "examples__${test_name}"
      "${CGAL_INSTALLATION_FILES_DIR}/examples"
      "${cgal_config_path}"
      TEST_EXAMPLES
      ${_recursive_call_extra_args}
    )
  endif()
endfunction()

##
# Function: test_CGAL_installation
# Adds config and install tests for a CGAL installation
##
function(test_CGAL_installation test_name install_prefix build_dir)
  add_test(NAME config_${test_name}
    COMMAND ${CMAKE_COMMAND} ${GENERATOR} ${INIT_FILE} "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
                                                       "-DCMAKE_INSTALL_PREFIX=${install_prefix}"
                                                       "-DCGAL_INSTALL_LIB_DIR=lib"
                                                       -S "${CGAL_SOURCE_DIR}" -B "${build_dir}")
  add_test(NAME install_${test_name}
    COMMAND ${CMAKE_COMMAND} --install "${build_dir}" --config "$<CONFIG>")
  set_property(TEST install_${test_name} APPEND PROPERTY DEPENDS config_${test_name})
  set_property(TEST install_${test_name}
                    config_${test_name}
               APPEND PROPERTY LABELS Installation_Tests CGAL_cmake_testsuite)
  set_property(TEST install_${test_name} APPEND PROPERTY FIXTURES_SETUP install_${test_name}_fixture)
endfunction()

CGAL_add_config_file_test(
  test_config_file_in_CGAL_SOURCE_DIR
  "${CMAKE_BINARY_DIR}/test_config_file"
  "${CGAL_SOURCE_DIR}/CGALConfig.cmake"
)

CGAL_add_config_file_test(
  test_config_file_in_lib_cmake_CGAL
  "${CMAKE_BINARY_DIR}/test_config_file"
  "${CGAL_INSTALLATION_FILES_DIR}/lib/cmake/CGAL/CGALConfig.cmake"
)

# install cgal in the non standard place without Qt6
get_filename_component(CORRECT_INSTALL_PATH "${WORKING_DIRS_PREFIX}/CGAL_cmake_test_suite/non_standard_install/dummy.txt" DIRECTORY)
test_CGAL_installation(
  non_standard_cgal
  "${CORRECT_INSTALL_PATH}"
  "${WORKING_DIRS_PREFIX}/CGAL_cmake_test_suite/non_standard_build"
)

#test CGAL_DIR=non standard place without cgal_qt
CGAL_add_config_file_test(
  test_config_file_non_standard_install_lib_cmake
  "${CMAKE_BINARY_DIR}/test_config_file"
  "${WORKING_DIRS_PREFIX}/CGAL_cmake_test_suite/non_standard_install/lib/cmake/CGAL/CGALConfig.cmake"
  FIXTURES_REQUIRED install_non_standard_cgal_fixture
)

if(CGAL_Qt6_FOUND)
  #configure cgal for a non standard install with Qt6
  test_CGAL_installation(
    non_standard_cgal_qt
    "${WORKING_DIRS_PREFIX}/CGAL_cmake_test_suite/non_standard_install_qt"
    "${WORKING_DIRS_PREFIX}/CGAL_cmake_test_suite/non_standard_build_qt"
  )

  #test CGAL_DIR=non standard place with cgal_qt
  CGAL_add_config_file_test(
    test_config_file_non_standard_install_lib_cmake_with_qt
    "${CMAKE_BINARY_DIR}/test_config_file_qt"
    "${WORKING_DIRS_PREFIX}/CGAL_cmake_test_suite/non_standard_install_qt/lib/cmake/CGAL/CGALConfig.cmake"
    FIXTURES_REQUIRED install_non_standard_cgal_qt_fixture
  )
endif()#CGAL_Qt6_FOUND

CGAL_add_config_file_test(
  test_config_file_non_standard
  "${CMAKE_BINARY_DIR}/test_config_file"
  "${WORKING_DIRS_PREFIX}/CGAL_cmake_test_suite/non_standard_build/CGALConfig.cmake"
  FIXTURES_REQUIRED install_non_standard_cgal_fixture
)
