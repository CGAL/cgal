# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Apollonius_graph_2_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( ag2_exact_traits.cpp )
  create_single_source_cgal_program( ag2_exact_traits_sqrt.cpp )
  create_single_source_cgal_program( ag2_filtered_traits_no_hidden.cpp )
  create_single_source_cgal_program( ag2_hierarchy.cpp )

endif()
