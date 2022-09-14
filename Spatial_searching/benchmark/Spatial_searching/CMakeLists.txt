# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Spatial_searching_)

find_package(CGAL REQUIRED COMPONENTS Core)

include_directories(BEFORE "include")

find_package(Eigen3 3.1.91 QUIET) # (requires 3.1.91 or greater)
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: These benchmarks require Eigen 3.1.91 (or greater), and will not be compiled.")
  return()
endif()

# create_single_source_cgal_program("Compare_ANN_STANN_CGAL.cpp") # does not compile, missing dependency
create_single_source_cgal_program("nanoflan.cpp")
create_single_source_cgal_program("binary.cpp")
create_single_source_cgal_program("nearest_neighbor_searching_50.cpp")
create_single_source_cgal_program("nearest_neighbor_searching_inplace_50.cpp")
create_single_source_cgal_program("Nearest_neighbor_searching.cpp")
create_single_source_cgal_program("Nearest_neighbor_searching_2D.cpp")
create_single_source_cgal_program("Nearest_neighbor_searching_2D_user_defined.cpp")
create_single_source_cgal_program("Split_data.cpp")
create_single_source_cgal_program("nn3cgal.cpp")
# create_single_source_cgal_program("nn4cgal.cpp") # this file does not exist for some reason
create_single_source_cgal_program("nn3nanoflan.cpp")
create_single_source_cgal_program("sizeof.cpp")
# create_single_source_cgal_program("deque.cpp") # does not compile, lots of errors
foreach(
  target
  # Compare_ANN_STANN_CGAL # see above
  nanoflan
  binary
  nearest_neighbor_searching_50
  nearest_neighbor_searching_inplace_50
  Nearest_neighbor_searching
  Nearest_neighbor_searching_2D
  Nearest_neighbor_searching_2D_user_defined
  Split_data
  nn3cgal
  # nn4cgal # see above
  nn3nanoflan
  sizeof
  # deque # see above
  )
  target_link_libraries(${target} PUBLIC CGAL::Eigen3_support)
endforeach()
