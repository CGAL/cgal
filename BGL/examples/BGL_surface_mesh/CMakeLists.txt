cmake_minimum_required(VERSION 3.1...3.15)
project(BGL_surface_mesh_Examples)

find_package(CGAL REQUIRED)

# include for local package

create_single_source_cgal_program("prim.cpp")
create_single_source_cgal_program("gwdwg.cpp")
create_single_source_cgal_program("seam_mesh.cpp")
create_single_source_cgal_program("write_inp.cpp")
create_single_source_cgal_program("surface_mesh_dual.cpp")
create_single_source_cgal_program("connected_components.cpp")

find_package(METIS)

if(METIS_FOUND)
  include_directories(${METIS_INCLUDE_DIRS})

  create_single_source_cgal_program("surface_mesh_partition.cpp")
  target_link_libraries(surface_mesh_partition PRIVATE ${METIS_LIBRARIES})
else()
  message(STATUS "Examples that use the METIS library will not be compiled.")
endif()
