# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Kernel_23_Test)

cmake_minimum_required(VERSION 2.6.2)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CGAL_CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)
  include_directories (BEFORE ../../../Cartesian_kernel/include ../../../Homogeneous_kernel/include ../../../Filtered_kernel/include)
  include_directories (BEFORE ../../include)

  create_single_source_cgal_program( Cartesian.cpp )
  create_single_source_cgal_program( Dimension.cpp )
  create_single_source_cgal_program( Filtered_cartesian.cpp )
  create_single_source_cgal_program( Filtered_homogeneous.cpp )
  create_single_source_cgal_program( Homogeneous.cpp )
  create_single_source_cgal_program( Kernel_checker.cpp )
  create_single_source_cgal_program( Lazy_kernel.cpp )
  create_single_source_cgal_program( Simple_cartesian.cpp )
  create_single_source_cgal_program( Simple_homogeneous.cpp )
  create_single_source_cgal_program( test_kernel__.cpp )

endif()
