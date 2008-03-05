# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Skin_surface_3_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( skin_surface_pdb_reader.cpp )
  create_single_source_cgal_program( skin_surface_simple.cpp )
  create_single_source_cgal_program( skin_surface_subdiv.cpp )
  create_single_source_cgal_program( skin_surface_subdiv_with_normals.cpp )
  create_single_source_cgal_program( union_of_balls_simple.cpp )
  create_single_source_cgal_program( union_of_balls_subdiv.cpp )

endif()
