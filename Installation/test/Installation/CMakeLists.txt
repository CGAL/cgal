# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.20)

project( Installation_Tests )


macro(create_link_to_program COMPONENT)

  add_executable(link_to_${COMPONENT} link_to_${COMPONENT}.cpp)

  include_directories(${${COMPONENT}_3RD_PARTY_INCLUDE_DIRS})

  add_definitions(${${COMPONENT}_3RD_PARTY_DEFINITIONS})

  link_directories(${${COMPONENT}_3RD_PARTY_LIBRARIES_DIRS})

  # Link the executable to CGAL and third-party libraries
  if(CGAL_AUTO_LINK_ENABLED OR CGAL_HEADER_ONLY)
    target_link_libraries(link_to_${COMPONENT} ${CGAL_3RD_PARTY_LIBRARIES}
                          ${${COMPONENT}_3RD_PARTY_LIBRARIES})
  else()
    target_link_libraries(
      link_to_${COMPONENT} CGAL::${COMPONENT} ${CGAL_3RD_PARTY_LIBRARIES}
      ${${COMPONENT}_3RD_PARTY_LIBRARIES})
  endif()

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS link_to_${COMPONENT})

endmacro()

find_package(CGAL REQUIRED COMPONENTS Core)

include(${CGAL_USE_FILE})

include(CGAL_CreateSingleSourceCGALProgram)

create_single_source_cgal_program("endian.cpp")

create_single_source_cgal_program("deprecation_warning.cpp")

create_single_source_cgal_program("test_use_h.cpp")

create_single_source_cgal_program(will_fail.cpp will_fail_aux.cpp)
if(BUILD_TESTING)
  set_property(TEST execution___of__will_fail PROPERTY WILL_FAIL TRUE)
endif()

find_package(TBB QUIET)
include(CGAL_TBB_support)
create_single_source_cgal_program("test_TBB.cpp")
if(TARGET CGAl::TBB_support)
  target_link_libraries(test_TBB PUBLIC CGAL::TBB_support)
endif()

create_link_to_program(CGAL)

if ( CGAL_Core_FOUND )
  create_link_to_program(CGAL_Core)
endif()

if(WITH_CGAL_ImageIO)
  create_link_to_program(CGAL_ImageIO)
endif()


if(WIN32 OR CMAKE_SYSTEM_NAME STREQUAL Windows)
  add_executable(display_dll_version_info display_dll_version_info.cpp)
  target_link_libraries(display_dll_version_info version)
  add_executable(test_gmp_mpfr_dll test_gmp_mpfr_dll.cpp)
  target_link_libraries(test_gmp_mpfr_dll version)
  CGAL_add_test(test_gmp_mpfr_dll)
  add_to_cached_list( CGAL_EXECUTABLE_TARGETS test_gmp_mpfr_dll )
  string(RANDOM RDM_DIR)#5 random chars to avoid conflicts in parallel testsuites
  set(NON_STANDARD_INSTALL_PREFIX $ENV{TEMP}/${RDM_DIR})
else()
  set(NON_STANDARD_INSTALL_PREFIX ${CMAKE_BINARY_DIR})
endif()

if(WITH_CGAL_Qt5)
  find_package(Qt5 QUIET COMPONENTS Core)
  if(Qt5_FOUND)
    create_link_to_program(CGAL_Qt5)
  endif()
endif()

find_package(LEDA QUIET)
if(LEDA_FOUND)
  include(${LEDA_USE_FILE})
  add_executable(test_LEDA_IO test_LEDA_IO.cpp)
  target_link_libraries(test_LEDA_IO ${LEDA_LIBRARIES})
endif()

if("$ENV{CMAKE_GENERATOR}" STREQUAL "")
    set(GENERATOR )
  else()
    set(GENERATOR $ENV{CMAKE_GENERATOR})
  endif()

  if("$ENV{INIT_FILE}" STREQUAL "")
  set(INIT_FILE )
  else()
  set(INIT_FILE -C$ENV{INIT_FILE})
  endif()

function(CGAL_installation_test_find_package_version mode)
  set(EXACT)
  if(mode STREQUAL "less")
    set(CGAL_MINOR_VERSION 0)
  endif()
  if(mode STREQUAL "less_major")
    set(CGAL_MAJOR_VERSION 0)
  endif()
  if(mode STREQUAL "greater" OR mode STREQUAL "fail-exact")
    math(EXPR CGAL_MINOR_VERSION "${CGAL_MINOR_VERSION} + 1")
  endif()
  if(mode STREQUAL "exact" OR mode STREQUAL "fail-exact")
    set(EXACT "EXACT ")
  endif()
  if(CGAL_BUGFIX_VERSION AND CGAL_BUGFIX_VERSION GREATER 0)
    set(VERSION
        "${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}.${CGAL_BUGFIX_VERSION}")
  else()
    set(VERSION "${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}")
  endif()
  file(REMOVE_RECURSE "${CMAKE_CURRENT_BINARY_DIR}/build-test_find_package_version_${mode}")
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build-test_find_package_version_${mode})
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_find_package_version_${mode})

  configure_file(test_find_package.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/test_find_package_version_${mode}/CMakeLists.txt)

  add_test(NAME test_find_package_version_${mode}
    COMMAND ${CMAKE_COMMAND} ${GENERATOR} ${INIT_FILE} -DCMAKE_BUILD_TYPE=Release ${CMAKE_CURRENT_BINARY_DIR}/test_find_package_version_${mode}
    )
  cgal_setup_test_properties(test_find_package_version_${mode})
  set_property(TEST test_find_package_version_${mode}
    PROPERTY
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build-test_find_package_version_${mode})
endfunction()

CGAL_installation_test_find_package_version(less)
CGAL_installation_test_find_package_version(less_major)
CGAL_installation_test_find_package_version(equal)
CGAL_installation_test_find_package_version(greater)
CGAL_installation_test_find_package_version(exact)
CGAL_installation_test_find_package_version(fail-exact)

set_tests_properties(
  test_find_package_version_greater
  test_find_package_version_less_major
  test_find_package_version_fail-exact
  PROPERTIES WILL_FAIL TRUE)

find_package(CGAL COMPONENTS Qt5)

file(MAKE_DIRECTORY "${NON_STANDARD_INSTALL_PREFIX}/non_standard_install")
file(MAKE_DIRECTORY ${NON_STANDARD_INSTALL_PREFIX}/non_standard_build)
file(MAKE_DIRECTORY ${NON_STANDARD_INSTALL_PREFIX}/test_config_file)

if(CGAL_Qt5_FOUND)
  file(MAKE_DIRECTORY ${NON_STANDARD_INSTALL_PREFIX}/non_standard_install_qt5)
  file(MAKE_DIRECTORY ${NON_STANDARD_INSTALL_PREFIX}/non_standard_build_qt5)
  file(MAKE_DIRECTORY ${NON_STANDARD_INSTALL_PREFIX}/test_config_file_qt5)
  configure_file(test_configuration_qt5.cmake.in ${CMAKE_BINARY_DIR}/test_config_file_qt5/CMakeLists.txt @ONLY)
endif()#CGAL_Qt5_FOUND

#If ctest is ran from a global config, CGAL_SOURCE_DIR exists, but from Installation/test it doesn't. In that case, however, there is a CGAL_DIR.
if("${CGAL_SOURCE_DIR}" STREQUAL "")
  if("$ENV{CGAL_RELEASE_DIR}" STREQUAL "")
    set(CGAL_SOURCE_DIR ${CGAL_DIR})
  else()
    set(CGAL_SOURCE_DIR $ENV{CGAL_RELEASE_DIR})
  endif()
endif()
string(REGEX REPLACE "/cygdrive/c/" "C:/" CGAL_SOURCE_DIR ${CGAL_SOURCE_DIR})

configure_file(test_configuration.cmake.in ${CMAKE_BINARY_DIR}/test_config_file/CMakeLists.txt @ONLY)

#test CGAL_DIR = source_dir (Git_root or CGAL-5.x dir.
get_filename_component(CGAL_DIR_CORRECT_PATH "${CGAL_SOURCE_DIR}/CGALConfig.cmake" DIRECTORY)
add_test(NAME test_config_file
  COMMAND ${CMAKE_COMMAND} ${GENERATOR} ${INIT_FILE} -DCMAKE_BUILD_TYPE=Release -DCGAL_DIR=${CGAL_DIR_CORRECT_PATH} -DCGAL_GIVEN_DIR=${CGAL_DIR_CORRECT_PATH}
    "${CMAKE_BINARY_DIR}/test_config_file" #src
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/build-test_config_file") #build
list(APPEND test_config_lst "test_config_file")

if(RUNNING_CGAL_AUTO_TEST OR CGAL_TEST_SUITE)#CGAL_BRANCH_BUILD
  #test CGAL_DIR = CGAL-5.x/lib/cmake/CGAL
  get_filename_component(CGAL_DIR_CORRECT_PATH "${CGAL_SOURCE_DIR}/lib/cmake/CGAL/CGALConfig.cmake" DIRECTORY)
  add_test(NAME test_config_file_2
    COMMAND ${CMAKE_COMMAND} ${GENERATOR} ${INIT_FILE} -DCMAKE_BUILD_TYPE=Release -DCGAL_DIR=${CGAL_DIR_CORRECT_PATH} -DCGAL_GIVEN_DIR=${CGAL_DIR_CORRECT_PATH}
      "${CMAKE_BINARY_DIR}/test_config_file"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/build-test_config_file_2")
else()#CGAL_BRANCH_BUILD
  #use the CGAL_DIR
  get_filename_component(CGAL_DIR_CORRECT_PATH "${CGAL_SOURCE_DIR}/Installation/lib/cmake/CGAL/CGALConfig.cmake" DIRECTORY)
  add_test(NAME test_config_file_2
    COMMAND ${CMAKE_COMMAND} ${GENERATOR} ${INIT_FILE} -DCMAKE_BUILD_TYPE=Release -DCGAL_DIR=${CGAL_DIR_CORRECT_PATH} -DCGAL_GIVEN_DIR=${CGAL_DIR_CORRECT_PATH}
      "${CMAKE_BINARY_DIR}/test_config_file"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/build-test_config_file_2")
endif()#CGAL_BRANCH_BUILD
list(APPEND test_config_lst "test_config_file_2")

#configure cgal for a non standard install without qt5
get_filename_component(CORRECT_INSTALL_PATH "${NON_STANDARD_INSTALL_PREFIX}/non_standard_install/dummy.txt" DIRECTORY)
add_test(NAME config_non_standard_cgal
  COMMAND ${CMAKE_COMMAND} ${GENERATOR} ${INIT_FILE} -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${CORRECT_INSTALL_PATH} -DCGAL_INSTALL_LIB_DIR=lib
  "${CGAL_SOURCE_DIR}"
  WORKING_DIRECTORY "${NON_STANDARD_INSTALL_PREFIX}/non_standard_build")
#install cgal in the non standard place
add_test(NAME install_non_standard_cgal
  COMMAND ${CMAKE_COMMAND} --build "${NON_STANDARD_INSTALL_PREFIX}/non_standard_build" --target "install" --config "$<CONFIG>")

#test CGAL_DIR=non standard place without cgal_qt5
get_filename_component(CGAL_DIR_CORRECT_PATH "${NON_STANDARD_INSTALL_PREFIX}/non_standard_install/lib/cmake/CGAL/CGALConfig.cmake" DIRECTORY)
add_test(NAME test_config_file_3
COMMAND ${CMAKE_COMMAND} ${GENERATOR} ${INIT_FILE} -DCMAKE_BUILD_TYPE=Release -DCGAL_DIR=${CGAL_DIR_CORRECT_PATH} -DCGAL_GIVEN_DIR=${CGAL_DIR_CORRECT_PATH}
"${CMAKE_BINARY_DIR}/test_config_file"
WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/build-test_config_file_3")
list(APPEND test_config_lst "test_config_file_3")

if(CGAL_Qt5_FOUND)
  #configure cgal for a non standard install with qt5
  add_test(NAME config_non_standard_cgal_qt5
    COMMAND ${CMAKE_COMMAND} ${GENERATOR} ${INIT_FILE} -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${NON_STANDARD_INSTALL_PREFIX}/non_standard_install_qt5 -DCGAL_INSTALL_LIB_DIR=lib
    "${CGAL_SOURCE_DIR}"
    WORKING_DIRECTORY "${NON_STANDARD_INSTALL_PREFIX}/non_standard_build_qt5")
  #install cgal in the non standard place
  add_test(NAME install_non_standard_cgal_qt5
    COMMAND ${CMAKE_COMMAND} --build "${NON_STANDARD_INSTALL_PREFIX}/non_standard_build_qt5" --target "install" --config "$<CONFIG>")

  #test CGAL_DIR=non standard place with cgal_qt5
  get_filename_component(CGAL_DIR_CORRECT_PATH "${NON_STANDARD_INSTALL_PREFIX}/non_standard_install_qt5/lib/cmake/CGAL/CGALConfig.cmake" DIRECTORY)
  add_test(NAME test_config_file_4
    COMMAND ${CMAKE_COMMAND} ${GENERATOR} ${INIT_FILE} -DCMAKE_BUILD_TYPE=Release -DCGAL_DIR=${CGAL_DIR_CORRECT_PATH} -DCGAL_GIVEN_DIR=${CGAL_DIR_CORRECT_PATH}
    "${CMAKE_BINARY_DIR}/test_config_file_qt5"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/build-test_config_file_4")
  list(APPEND test_config_lst "test_config_file_4")
endif()#CGAL_Qt5_FOUND

#test CGAL_DIR=non standard build
get_filename_component(CGAL_DIR_CORRECT_PATH "${NON_STANDARD_INSTALL_PREFIX}/non_standard_build/CGALConfig.cmake" DIRECTORY)
add_test(NAME test_config_file_5
COMMAND ${CMAKE_COMMAND} ${GENERATOR} ${INIT_FILE} -DCMAKE_BUILD_TYPE=Release -DCGAL_DIR=${CGAL_DIR_CORRECT_PATH} -DCGAL_GIVEN_DIR=${CGAL_DIR_CORRECT_PATH}
 "${CMAKE_BINARY_DIR}/test_config_file"
 WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/build-test_config_file_5")
list(APPEND test_config_lst "test_config_file_5")

foreach(tgt ${test_config_lst})
  #add_custom_target(${tgt}_target)
  add_test(NAME "${tgt}_SetupFixture"  COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/build-${tgt}")
  add_test(NAME "${tgt}_CleanupFixture"  COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}/build-${tgt}")
  add_test(NAME compile__${tgt}
    COMMAND ${CMAKE_COMMAND} --build "${CMAKE_BINARY_DIR}/build-${tgt}" --config "$<CONFIG>")
    add_test(NAME exec__${tgt}
      COMMAND ./test_configuration
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/build-${tgt}")

  set_property(TEST "${tgt}_SetupFixture" APPEND PROPERTY FIXTURES_SETUP ${tgt}_target)
  set_property(TEST "${tgt}_CleanupFixture" APPEND PROPERTY FIXTURES_CLEANUP ${tgt}_target)
  set_property(TEST ${tgt} APPEND PROPERTY FIXTURES_REQUIRED ${tgt}_target)
  set_property(TEST compile__${tgt} APPEND PROPERTY FIXTURES_REQUIRED ${tgt}_target)
  set_property(TEST compile__${tgt} APPEND PROPERTY DEPENDS ${tgt})
  set_property(TEST exec__${tgt} APPEND PROPERTY FIXTURES_REQUIRED ${tgt}_target)
  set_property(TEST exec__${tgt} APPEND PROPERTY DEPENDS compile__${tgt})

  set_property(TEST "${tgt}_SetupFixture" APPEND PROPERTY LABELS Installation_Tests CGAL_cmake_testsuite)
  set_property(TEST "${tgt}_CleanupFixture" APPEND PROPERTY LABELS Installation_Tests CGAL_cmake_testsuite)
  set_property(TEST ${tgt} APPEND PROPERTY LABELS Installation_Tests CGAL_cmake_testsuite)
  set_property(TEST compile__${tgt} APPEND PROPERTY LABELS Installation_Tests CGAL_cmake_testsuite)
  set_property(TEST exec__${tgt} APPEND PROPERTY LABELS Installation_Tests CGAL_cmake_testsuite)

endforeach()

set_property(TEST install_non_standard_cgal APPEND PROPERTY DEPENDS config_non_standard_cgal)
set_property(TEST config_non_standard_cgal APPEND PROPERTY FIXTURES_SETUP test_config_file_3_target)
set_property(TEST install_non_standard_cgal APPEND PROPERTY FIXTURES_SETUP test_config_file_3_target)
set_property(TEST test_config_file_3 test_config_file_5 APPEND PROPERTY FIXTURES_REQUIRED test_config_file_3_target)

set_property(TEST install_non_standard_cgal APPEND PROPERTY LABELS Installation_Tests CGAL_cmake_testsuite)
set_property(TEST config_non_standard_cgal APPEND PROPERTY LABELS Installation_Tests CGAL_cmake_testsuite)

if(CGAL_Qt5_FOUND)
  set_property(TEST install_non_standard_cgal_qt5 APPEND PROPERTY DEPENDS config_non_standard_cgal_qt5)
  set_property(TEST config_non_standard_cgal_qt5 APPEND PROPERTY FIXTURES_SETUP test_config_file_4_target)
  set_property(TEST install_non_standard_cgal_qt5 APPEND PROPERTY FIXTURES_SETUP test_config_file_4_target)
  set_property(TEST test_config_file_4 APPEND PROPERTY FIXTURES_REQUIRED test_config_file_4_target)

  set_property(TEST install_non_standard_cgal_qt5 APPEND PROPERTY LABELS Installation_Tests CGAL_cmake_testsuite)
  set_property(TEST config_non_standard_cgal_qt5 APPEND PROPERTY LABELS Installation_Tests CGAL_cmake_testsuite)
endif()#CGAL_Qt5_FOUND
