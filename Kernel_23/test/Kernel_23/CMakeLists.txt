# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Kernel_23_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( Cartesian.cpp )
  create_single_source_cgal_program( Dimension.cpp )
  create_single_source_cgal_program( Filtered_cartesian.cpp )
  create_single_source_cgal_program( Filtered_homogeneous.cpp )
  create_single_source_cgal_program( Homogeneous.cpp )
  create_single_source_cgal_program( Kernel_checker.cpp )
  create_single_source_cgal_program( Lazy_kernel.cpp )
  create_single_source_cgal_program( Simple_cartesian.cpp )
  create_single_source_cgal_program( Simple_homogeneous.cpp )
  create_single_source_cgal_program( test_kernel_archetype_2.cpp )
  create_single_source_cgal_program( test_kernel_archetype_3.cpp )
  create_single_source_cgal_program( test_kernel__.cpp )

endif()
