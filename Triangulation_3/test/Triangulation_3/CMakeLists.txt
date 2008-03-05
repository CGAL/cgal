# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Triangulation_3_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( test_delaunay_3.cpp )
  create_single_source_cgal_program( test_delaunay_3_old_remove.cpp )
  create_single_source_cgal_program( test_delaunay_hierarchy_3.cpp )
  create_single_source_cgal_program( test_regular_3.cpp )
  create_single_source_cgal_program( test_regular_as_delaunay_3.cpp )
  create_single_source_cgal_program( test_regular_remove_3.cpp )
  create_single_source_cgal_program( test_regular_traits_3.cpp )
  create_single_source_cgal_program( test_simplex_3.cpp )
  create_single_source_cgal_program( test_triangulation_3.cpp )
  create_single_source_cgal_program( test_triangulation_tds.cpp )

endif()
