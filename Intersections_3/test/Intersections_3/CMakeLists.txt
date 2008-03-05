# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Intersections_3_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( isocub_isocub.cpp )
  create_single_source_cgal_program( line_isocub.cpp )
  create_single_source_cgal_program( plane_line.cpp )
  create_single_source_cgal_program( plane_plane.cpp )
  create_single_source_cgal_program( plane_ray.cpp )
  create_single_source_cgal_program( plane_seg.cpp )
  create_single_source_cgal_program( ray_isocub.cpp )
  create_single_source_cgal_program( seg_isocub.cpp )
  create_single_source_cgal_program( triangle_other.cpp )

endif()
