# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Kernel_d_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( afftrafo-test.cpp )
  create_single_source_cgal_program( interface-test.cpp )
  create_single_source_cgal_program( intersection-test.cpp )
  create_single_source_cgal_program( Iso_box_d_with_Cartesian.cpp )
  create_single_source_cgal_program( Iso_box_d_with_homogeneous.cpp )
  create_single_source_cgal_program( Linear_algebra-test.cpp )

endif()
