# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Surface_mesher_Examples)

find_package(CGAL REQUIRED COMPONENTS ImageIO)

if(CGAL_ImageIO_FOUND)

  create_single_source_cgal_program("mesh_a_3d_gray_image.cpp")
  create_single_source_cgal_program("mesh_an_implicit_function.cpp")

else()
  if(RUNNING_CGAL_AUTO_TEST OR CGAL_TEST_SUITE)
    # Just to avoid a warning from CMake if that variable is set on the command line...
  endif()

  message("NOTICE: This project requires the CGAL and CGAL ImageIO libraries, and will not be compiled.")

endif()
