# This is the CMake script for compiling this folder.

cmake_minimum_required(VERSION 3.1...3.23)
project(Poisson_surface_reconstruction_3)

if(CGAL_ENABLE_TESTING)
  enable_testing()
endif()

add_subdirectory("examples/Poisson_surface_reconstruction_3")
add_subdirectory("test/Poisson_surface_reconstruction_3")
