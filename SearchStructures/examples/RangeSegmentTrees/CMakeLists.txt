# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (RangeSegmentTrees_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( range_tree_1.cpp )
  create_single_source_cgal_program( range_tree_2.cpp )
  create_single_source_cgal_program( range_tree_3.cpp )
  create_single_source_cgal_program( range_tree_4.cpp )
  create_single_source_cgal_program( range_tree_map_2.cpp )
  create_single_source_cgal_program( range_tree_set_2.cpp )
  create_single_source_cgal_program( segment_tree_1.cpp )
  create_single_source_cgal_program( segment_tree_2.cpp )
  create_single_source_cgal_program( segment_tree_3.cpp )
  create_single_source_cgal_program( segment_tree_4.cpp )
  create_single_source_cgal_program( segment_tree_map_2.cpp )
  create_single_source_cgal_program( segment_tree_set_2.cpp )
  create_single_source_cgal_program( segment_tree_set_3.cpp )

endif()
