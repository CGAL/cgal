# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGALROOT}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( pdb_cat.cpp )
  create_single_source_cgal_program( pdb_echo.cpp )
  create_single_source_cgal_program( pdb_split.cpp )
  create_single_source_cgal_program( pdb_transform.cpp )
  create_single_source_cgal_program( pdb_tree.cpp )

endif()
