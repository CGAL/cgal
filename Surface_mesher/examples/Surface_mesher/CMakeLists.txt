# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Surface_mesher_Examples)

find_package(CGAL REQUIRED COMPONENTS ImageIO)

if(CGAL_ImageIO_FOUND)

  create_single_source_cgal_program("mesh_a_3d_gray_image.cpp")
  create_single_source_cgal_program("mesh_an_implicit_function.cpp")

else()
  if(RUNNING_CGAL_AUTO_TEST)
    # Just to avoid a warning from CMake if that variable is set on the command line...
  endif()

  message(
    STATUS
      "NOTICE: This program requires the CGAL and CGAL ImageIO libraries, and will not be compiled."
  )

endif()
