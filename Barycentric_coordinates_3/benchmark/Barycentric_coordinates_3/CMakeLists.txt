cmake_minimum_required(VERSION 3.12...3.31)
project(Barycentric_coordinates_3_Examples)

set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} -Wall -Wextra -Werror")

set(CMAKE_CXX_STANDARD 14)

find_package(CGAL QUIET COMPONENTS Core)
if(CGAL_FOUND)

  include(${CGAL_USE_FILE})
  include(CGAL_CreateSingleSourceCGALProgram)

  create_single_source_cgal_program("benchmark_tetrahedon_coordinates.cpp")
  create_single_source_cgal_program("benchmark_polyhedron_8_vertices.cpp")

else()
  message(WARNING "This program requires the CGAL library, and will not be compiled.")
endif()
