# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( comparison.cpp )
  create_single_source_cgal_program( complex_construction.cpp )
  create_single_source_cgal_program( construction.cpp )
  create_single_source_cgal_program( exploration_SM.cpp )
  create_single_source_cgal_program( extended_kernel.cpp )
  create_single_source_cgal_program( handling_double_coordinates.cpp )
  create_single_source_cgal_program( interface_polyhedron.cpp )
  create_single_source_cgal_program( nefIO.cpp )
  create_single_source_cgal_program( offIO.cpp )
  create_single_source_cgal_program( point_location.cpp )
  create_single_source_cgal_program( point_set_operations.cpp )
  create_single_source_cgal_program( polygon_construction.cpp )
  create_single_source_cgal_program( polyline_construction.cpp )
  create_single_source_cgal_program( set_operations.cpp )
  create_single_source_cgal_program( shell_exploration.cpp )
  create_single_source_cgal_program( simple.cpp )
  create_single_source_cgal_program( topological_operations.cpp )
  create_single_source_cgal_program( transformation.cpp )

endif()
