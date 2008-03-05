# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Arrangement_2_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( construction_test_suite_generator.cpp )
  create_single_source_cgal_program( ex_kernel_point.cpp )
  create_single_source_cgal_program( ex_kernel_segment.cpp )
  create_single_source_cgal_program( point_location.cpp )
  create_single_source_cgal_program( test_compute_curve_zone.cpp )
  create_single_source_cgal_program( test_construction.cpp )
  create_single_source_cgal_program( test_do_intersect.cpp )
  create_single_source_cgal_program( test_dual.cpp )
  create_single_source_cgal_program( test_insertion.cpp )
  create_single_source_cgal_program( test_iso_verts.cpp )
  create_single_source_cgal_program( test_removal.cpp )
  create_single_source_cgal_program( test_traits.cpp )
  create_single_source_cgal_program( test_vert_ray_shoot_vert_segments.cpp )

endif()
