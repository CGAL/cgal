# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Straight_skeleton_2_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( Create_offset_polygons_2.cpp )
  create_single_source_cgal_program( Create_skeleton_and_offset_polygons_2.cpp )
  create_single_source_cgal_program( Create_skeleton_and_offset_polygons_from_polygon_with_holes_2.cpp )
  create_single_source_cgal_program( Create_skeleton_and_offset_polygons_with_holes_2.cpp )
  create_single_source_cgal_program( Create_straight_skeleton_2.cpp )
  create_single_source_cgal_program( Create_straight_skeleton_from_polygon_with_holes_2.cpp )
  create_single_source_cgal_program( Low_level_API.cpp )

endif()
