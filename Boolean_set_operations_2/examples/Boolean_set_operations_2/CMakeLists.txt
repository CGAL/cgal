# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( bezier_traits_adapter.cpp )
  create_single_source_cgal_program( circle_segment.cpp )
  create_single_source_cgal_program( conic_traits_adapter.cpp )
  create_single_source_cgal_program( connect_polygon.cpp )
  create_single_source_cgal_program( do_intersect.cpp )
  create_single_source_cgal_program( dxf_union.cpp )
  create_single_source_cgal_program( sequence.cpp )
  create_single_source_cgal_program( set_union.cpp )
  create_single_source_cgal_program( simple_join_intersect.cpp )
  create_single_source_cgal_program( symmetric_difference.cpp )

endif()
