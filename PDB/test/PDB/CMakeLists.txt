# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (PDB_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( check_hetatom.cpp )
  create_single_source_cgal_program( check_iterators.cpp )
  create_single_source_cgal_program( check_pdb.cpp )
  create_single_source_cgal_program( check_protein.cpp )
  create_single_source_cgal_program( check_small_map.cpp )

endif()
