# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Segment_Delaunay_graph_2_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( test_sdg_2_et.cpp )
  create_single_source_cgal_program( test_sdg_2_ft.cpp )
  create_single_source_cgal_program( test_sdg_hierarchy_2_et.cpp )
  create_single_source_cgal_program( test_sdg_hierarchy_2_ft.cpp )
  create_single_source_cgal_program( test_sdg_hierarchy_info_2_no_up.cpp )
  create_single_source_cgal_program( test_sdg_hierarchy_info_2_up.cpp )
  create_single_source_cgal_program( test_sdg_info_2.cpp )
  create_single_source_cgal_program( test_sdg_traits_2.cpp )

endif()
