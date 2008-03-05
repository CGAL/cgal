# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Nef_2_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( EPoint-test.cpp )
  create_single_source_cgal_program( Nef_polyhedron_2-test.cpp )
  create_single_source_cgal_program( point_location.cpp )
  create_single_source_cgal_program( Polynomial-test.cpp )

endif()
