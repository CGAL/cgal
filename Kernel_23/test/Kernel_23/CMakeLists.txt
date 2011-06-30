# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Kernel_23_Test)

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

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
  create_single_source_cgal_program( overload_bug.cpp )
endif()
