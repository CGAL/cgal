# Top level CMakeLists.txt for CGAL-branchbuild

# Minimal version of CMake:
cmake_minimum_required(VERSION 3.12...3.29)

message("== CMake setup ==")
project(CGAL CXX C)
export(PACKAGE CGAL)

set(CGAL_BRANCH_BUILD
    ON
    CACHE INTERNAL "Create CGAL from a Git branch" FORCE)

include(${CMAKE_CURRENT_SOURCE_DIR}/CGALConfigVersion.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/Installation/cmake/modules/CGAL_SCM.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/Installation/cmake/modules/CGAL_enable_end_of_configuration_hook.cmake)
cgal_detect_git(${CMAKE_CURRENT_SOURCE_DIR})

function(CGAL_error_if_detect_in_source_build)
  # If in a Git repository, forbid in-source builds
  get_filename_component(srcdir "${CMAKE_CURRENT_SOURCE_DIR}" REALPATH)
  get_filename_component(bindir "${CMAKE_CURRENT_BINARY_DIR}" REALPATH)
  if("${srcdir}" STREQUAL "${bindir}")
    message(
      FATAL_ERROR
        [=[
############
Since CGAL-4.12.1, you can no longer configure an in-source build in a Git
repository. See this StackOverlow question and answers for a way to create
a separate build directory:
  https://stackoverflow.com/q/45518317/1728537
############
]=])
  endif()
endfunction()

if("${CGAL_SCM_NAME}" STREQUAL "git")
  cgal_error_if_detect_in_source_build()
endif()

# add option for duplicate file detection
option(
  CGAL_REPORT_DUPLICATE_FILES
  "Switch on to start (naive) detection of duplicate source- and headerfiles in packages"
  OFF)

message("== CMake setup (DONE) ==\n")

# Enable testing with CGAL_ENABLE_TESTING. Before CGAL-6.0, users would enable
# the tests by specifying BUILD_TESTING. For compatibility, If BUILD_TESTING is
# set, that is the default value for CGAL_ENABLE_TESTING. Otherwise, the default
# value is OFF.
option(CGAL_ENABLE_TESTING "Build the testing tree." ${BUILD_TESTING})
if(CGAL_ENABLE_TESTING AND NOT POLICY CMP0064)
  message(
    FATAL_ERROR "CGAL support of CTest requires CMake version 3.4 or later.
The variable CGAL_ENABLE_TESTING must be set of OFF.")
endif()
if(CGAL_ENABLE_TESTING)
  enable_testing()
endif()

#setup prefix path
file(GLOB BOOST_DIRS "C:/Program Files/boost/*")
if(NOT BOOST_DIRS)
  file(GLOB BOOST_DIRS "C:/Program Files/boost*")
endif()
if(NOT BOOST_DIRS)
  file(GLOB BOOST_DIRS "C:/Program Files (x86)/boost/*")
endif()
if(NOT BOOST_DIRS)
  file(GLOB BOOST_DIRS "C:/Program Files (x86)/boost*")
endif()
if(NOT BOOST_DIRS)
  file(GLOB BOOST_DIRS "C:/local/boost*")
endif()
if(BOOST_DIRS)
  list(GET BOOST_DIRS 0 boost_dir)
  list(APPEND CMAKE_PREFIX_PATH "${boost_dir}")
endif()#BOOST_DIRS

file(GLOB QT_DIRS "C:/Qt/5")
if(QT_DIRS)
  list(GET QT_DIRS 0 qt_dir)
  file(GLOB COMPS "${qt_dir}/msvc*")
  if(COMPS)
    list(GET COMPS 0 COMP)
    list(APPEND CMAKE_PREFIX_PATH "${COMP}")
  endif()#COMPS
endif()#QT_DIRS


# and finally start actual build
add_subdirectory(Installation)
add_subdirectory(Documentation/doc)

if(NOT TARGET uninstall)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake IMMEDIATE @ONLY)

  add_custom_target(
    uninstall COMMAND ${CMAKE_COMMAND} -P
                      ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
