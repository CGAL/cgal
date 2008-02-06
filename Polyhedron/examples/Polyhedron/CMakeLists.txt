# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGALROOT}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( polyhedron_prog_color.cpp )
  create_single_source_cgal_program( polyhedron_prog_cube.cpp )
  create_single_source_cgal_program( polyhedron_prog_cut_cube.cpp )
  create_single_source_cgal_program( polyhedron_prog_incr_builder.cpp )
  create_single_source_cgal_program( polyhedron_prog_normals.cpp )
  create_single_source_cgal_program( polyhedron_prog_off.cpp )
  create_single_source_cgal_program( polyhedron_prog_planes.cpp )
  create_single_source_cgal_program( polyhedron_prog_simple.cpp )
  create_single_source_cgal_program( polyhedron_prog_subdiv.cpp )
  create_single_source_cgal_program( polyhedron_prog_subdiv_with_boundary.cpp )
  create_single_source_cgal_program( polyhedron_prog_tetra.cpp )
  create_single_source_cgal_program( polyhedron_prog_traits.cpp )
  create_single_source_cgal_program( polyhedron_prog_vector.cpp )
  create_single_source_cgal_program( polyhedron_prog_vertex_color.cpp )
  create_single_source_cgal_program( polyhedron_prog_vertex_normal.cpp )
  create_single_source_cgal_program( polyhedron_self_intersection.cpp )

endif()
