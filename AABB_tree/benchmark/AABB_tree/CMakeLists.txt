# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project( AABB_traits_benchmark)

# CGAL and its components
find_package( CGAL REQUIRED )
include( ${CGAL_USE_FILE} )

# Boost and its components
find_package( Boost REQUIRED )
# include for local directory

# google benchmark
find_package(benchmark REQUIRED)

# include for local package
include_directories( BEFORE ../../include )

# add_executable (test_ test.cpp) // TODO: fix this benchmark
add_executable(tree_creation tree_creation.cpp)
target_link_libraries(tree_creation benchmark::benchmark)

