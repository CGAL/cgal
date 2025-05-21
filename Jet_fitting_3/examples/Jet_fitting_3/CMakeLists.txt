# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Jet_fitting_3_Examples)

find_package(CGAL REQUIRED)

# use Eigen
find_package(Eigen3 3.1.0 QUIET) #(requires 3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)

  # Link with Boost.ProgramOptions (optional)
  find_package(Boost QUIET COMPONENTS program_options)
  if(Boost_PROGRAM_OPTIONS_FOUND)
    create_single_source_cgal_program("Mesh_estimation.cpp")
    target_link_libraries(Mesh_estimation PRIVATE CGAL::Eigen3_support
                                                  Boost::program_options)
  else()
    message(STATUS "NOTICE: The example 'Mesh_estimation' requires Boost Program Options, and will not be compiled.")
  endif()

  create_single_source_cgal_program("Single_estimation.cpp")
  target_link_libraries(Single_estimation PRIVATE CGAL::Eigen3_support)

else()
  message("NOTICE: These examples require Eigen 3.1 (or greater), and will not be compiled.")
endif()
