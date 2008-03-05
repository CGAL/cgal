# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Circular_kernel_3_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( test_Exact_spherical_kernel.cpp )
  create_single_source_cgal_program( test_Lazy_Spherical_kernel.cpp )
  create_single_source_cgal_program( test_Spherical_kernel.cpp )

endif()
