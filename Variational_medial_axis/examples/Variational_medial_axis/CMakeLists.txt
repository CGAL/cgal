# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project( Variational_medial_axis )

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core )
find_package(Eigen3 3.2.0 QUIET) #(requires 3.2.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
	message(STATUS "Found Eigen")
	create_single_source_cgal_program( "vmas_example.cpp" )
	target_link_libraries(vmas_example PRIVATE CGAL::Eigen3_support)
else()
  message(ERROR "This program requires the Eigen library, and will not be compiled.")
endif()


