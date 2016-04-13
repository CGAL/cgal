# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Circular_kernel_3_test )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

include_directories (BEFORE include)
include_directories (BEFORE ../../include)
include_directories (BEFORE ../../../Kernel_23/test/Kernel_23/include)
include_directories (BEFORE ../Kernel_23/include)
include_directories (BEFORE ../Cartesian_kernel/include)
include_directories (BEFORE ../Homogeneous_kernel/include)
include_directories (BEFORE ../../../Intersections_3/include)
include_directories (BEFORE ../../../Algebraic_kernel_for_spheres/include)
include_directories (BEFORE ../../../Kernel_23/include)
include_directories (BEFORE ../../../Cartesian_kernel/include)

create_single_source_cgal_program( "test_Spherical_kernel.cpp" )
create_single_source_cgal_program( "test_Spherical_kernel_basics.cpp" )
create_single_source_cgal_program( "test_Exact_spherical_kernel.cpp" )
create_single_source_cgal_program( "test_Lazy_Spherical_kernel.cpp" )
create_single_source_cgal_program( "test_Lazy_spherical_kernel_basics.cpp" )
create_single_source_cgal_program( "test_Spherical_kernel_with_core.cpp" )

else()
  
    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
  
endif()

