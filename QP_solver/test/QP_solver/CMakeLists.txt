# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGALROOT}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( master_mps_to_derivatives.cpp )
  create_single_source_cgal_program( test_default_bounds.cpp )
  create_single_source_cgal_program( test_MPS.cpp )
  create_single_source_cgal_program( test_random_qp2.cpp )
  create_single_source_cgal_program( test_random_qp.cpp )
  create_single_source_cgal_program( test_solver.cpp )

endif()
