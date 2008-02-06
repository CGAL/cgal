# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGALROOT}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( NGHK_skin_surface_simple.cpp )
  create_single_source_cgal_program( NGHK_skin_surface_subdiv.cpp )
  create_single_source_cgal_program( skin_surface_pdb_reader.cpp )
  create_single_source_cgal_program( skin_surface_simple.cpp )
  create_single_source_cgal_program( skin_surface_subdiv.cpp )
  create_single_source_cgal_program( skin_surface_subdiv_with_normals.cpp )
  create_single_source_cgal_program( skin_surface_with_surface_mesher.cpp )
  create_single_source_cgal_program( union_of_balls_simple.cpp )
  create_single_source_cgal_program( union_of_balls_subdiv.cpp )

endif()
