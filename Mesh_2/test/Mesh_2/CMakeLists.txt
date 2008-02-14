# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Mesh_2_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( bench_double_map.cpp )
  create_single_source_cgal_program( test_conforming.cpp )
  create_single_source_cgal_program( test_double_map.cpp )
  create_single_source_cgal_program( test_filtred_container.cpp )
  create_single_source_cgal_program( test_meshing.cpp )
  create_single_source_cgal_program( test_vertex_and_face_bases.cpp )

endif()
