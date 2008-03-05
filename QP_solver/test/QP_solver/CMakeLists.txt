# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (QP_solver_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( master_mps_to_derivatives.cpp )
  create_single_source_cgal_program( test_default_bounds.cpp )
  create_single_source_cgal_program( test_MPS.cpp )
  create_single_source_cgal_program( test_random_qp2.cpp )
  create_single_source_cgal_program( test_random_qp.cpp )
  create_single_source_cgal_program( test_solver.cpp )

endif()
