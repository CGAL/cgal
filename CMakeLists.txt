# Top level CMakeLists.txt for CGAL-branchbuild
message( "== CMake setup ==" )
project(CGAL CXX C)

# Minimal version of CMake:
cmake_minimum_required(VERSION 2.8.11)

set( CGAL_BRANCH_BUILD ON CACHE INTERNAL "Create CGAL from a Git branch" FORCE)

include(${CMAKE_SOURCE_DIR}/Installation/cmake/modules/CGAL_SCM.cmake)
CGAL_detect_git(${CMAKE_SOURCE_DIR})

# add option for duplicate file detection
option( CGAL_REPORT_DUPLICATE_FILES "Switch on to start (naive) detection of duplicate source- and headerfiles in packages" OFF)

message( "== CMake setup (DONE) ==\n" )

# Enable testing with BUILD_TESTING
option(BUILD_TESTING "Build the testing tree." OFF)
include(CTest)
if(BUILD_TESTING AND NOT POLICY CMP0064)
  message(FATAL_ERROR "CGAL support of CTest requires CMake version 3.4 or later.
The variable BUILD_TESTING must be set of OFF.")
endif()

# and finally start actual build
add_subdirectory( Installation )
add_subdirectory( Documentation/doc)

if(NOT TARGET uninstall)
  configure_file(
      ${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in
      ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
      IMMEDIATE @ONLY)

  add_custom_target(uninstall
      COMMAND ${CMAKE_COMMAND} -P
      ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
