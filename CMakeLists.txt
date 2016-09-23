# Top level CMakeLists.txt for CGAL-branchbuild
message( "== CMake setup ==" )
project(CGAL CXX C)

# Minimal version of CMake:
cmake_minimum_required(VERSION 2.8.11)

set( CGAL_BRANCH_BUILD ON CACHE INTERNAL "Create CGAL from a Git branch" FORCE)

include(${CMAKE_SOURCE_DIR}/Installation/cmake/modules/CGAL_SCM.cmake)

# add option for duplicate file detection
option( CGAL_REPORT_DUPLICATE_FILES "Switch on to start (naive) detection of duplicate source- and headerfiles in packages" OFF)

message( "== CMake setup (DONE) ==\n" )

# Enable testing with BUILD_TESTING
option(BUILD_TESTING "Build the testing tree." OFF)
include(CTest)

# and finally start actual build
add_subdirectory( Installation )
add_subdirectory( Documentation/doc)
