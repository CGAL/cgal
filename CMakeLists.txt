# Top level CMakeLists.txt for CGAL-branchbuild
message( "== CMake setup ==" )
project(CGAL CXX C)

# Minimal version of CMake:
cmake_minimum_required(VERSION 2.8.11)

# option for branch build

option( CGAL_BRANCH_BUILD "Create CGAL from branch" ON)
mark_as_advanced( CGAL_BRANCH_BUILD )

include(${CMAKE_SOURCE_DIR}/Installation/cmake/modules/CGAL_SCM.cmake)

# add option for duplicate file detection
option( CGAL_REPORT_DUPLICATE_FILES "Switch on to start (naive) detection of duplicate source- and headerfiles in packages" OFF)

message( "== CMake setup (DONE) ==\n" )
# and finally start actual build
add_subdirectory( Installation )
add_subdirectory( Documentation/doc)
