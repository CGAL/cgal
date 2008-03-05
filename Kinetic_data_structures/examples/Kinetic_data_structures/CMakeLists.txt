# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Kinetic_data_structures_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( Delaunay_triangulation_2.cpp )
  create_single_source_cgal_program( Delaunay_triangulation_3.cpp )
  create_single_source_cgal_program( insert_event.cpp )
  create_single_source_cgal_program( regular_triangulation_3.cpp )
  create_single_source_cgal_program( sort.cpp )
  create_single_source_cgal_program( sweepline.cpp )

endif()
