# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGALROOT}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( test_andreas_example.cpp )
  create_single_source_cgal_program( test_PD.cpp )
  create_single_source_cgal_program( test_Polytope_distance_d_2.cpp )
  create_single_source_cgal_program( test_Polytope_distance_d_3.cpp )
  create_single_source_cgal_program( test_Polytope_distance_d_d.cpp )
  create_single_source_cgal_program( wilms_bug.cpp )
  create_single_source_cgal_program( zwick_bug.cpp )

endif()
