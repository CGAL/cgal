# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Surface_mesher_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( combined_spheres.cpp )
  create_single_source_cgal_program( implicit_surface_mesher_test.cpp )
  if(EXISTS internal-test-surface-mesher.cpp)
    create_single_source_cgal_program( internal-test-surface-mesher.cpp )
  endif()
  create_single_source_cgal_program( test_c2t3_iterators.cpp )
  create_single_source_cgal_program( test_robust_circumcenter.cpp )

endif()
