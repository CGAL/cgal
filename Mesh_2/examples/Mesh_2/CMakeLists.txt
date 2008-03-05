# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Mesh_2_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( conforming.cpp )
  create_single_source_cgal_program( mesh_class.cpp )
  create_single_source_cgal_program( mesh_global.cpp )
  create_single_source_cgal_program( mesh_with_seeds.cpp )

endif()
