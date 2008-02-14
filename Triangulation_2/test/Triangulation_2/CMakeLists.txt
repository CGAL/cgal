# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Triangulation_2_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( test_const_del_triangulation_2.cpp )
  create_single_source_cgal_program( test_constrained_triangulation_2.cpp )
  create_single_source_cgal_program( test_const_triang_plus_2.cpp )
  create_single_source_cgal_program( test_delaunay_hierarchy_2.cpp )
  create_single_source_cgal_program( test_delaunay_triangulation_2.cpp )
  create_single_source_cgal_program( test_regular_hierarchy_2.cpp )
  create_single_source_cgal_program( test_regular_triangulation_2.cpp )
  create_single_source_cgal_program( test_triangulation_2_bis.cpp )
  create_single_source_cgal_program( test_triangulation_2.cpp )
  create_single_source_cgal_program( test_triangulation_geom_traits.cpp )
  create_single_source_cgal_program( test_triangulation_tds.cpp )

endif()
