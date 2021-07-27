# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

project(Barycentric_coordinates_3_Examples)
set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} -Wall -Wextra -Werror")

cmake_minimum_required(VERSION 3.1...3.15)
set(CMAKE_CXX_STANDARD 14)

find_package(CGAL QUIET COMPONENTS Core)
if(CGAL_FOUND)

  include(${CGAL_USE_FILE})
  include(CGAL_CreateSingleSourceCGALProgram)

  create_single_source_cgal_program("tetrahedron_coordinates.cpp")
  create_single_source_cgal_program("shape_deformation.cpp")

else()
  message(WARNING "This program requires the CGAL library, and will not be compiled.")
endif()
