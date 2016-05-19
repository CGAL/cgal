# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Circular_kernel_2_test )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

include_directories (BEFORE include)
include_directories (BEFORE ../../include)
include_directories (BEFORE ../../../Kernel_23/test/Kernel_23/include)
include_directories (BEFORE ../Kernel_23/include)
include_directories (BEFORE ../Cartesian_kernel/include)
include_directories (BEFORE ../Homogeneous_kernel/include)
include_directories (BEFORE ../Intersections_2/include)

create_single_source_cgal_program( "test_Circular_kernel.cpp" )
create_single_source_cgal_program( "test_Lazy_circular_kernel.cpp" )
create_single_source_cgal_program( "test_Exact_circular_kernel.cpp" )
create_single_source_cgal_program( "test_Filtered_bbox_circular_kernel.cpp" )
create_single_source_cgal_program( "test_Line_arc.cpp" )
create_single_source_cgal_program( "test_Circular_kernel_basic.cpp" )
create_single_source_cgal_program( "test_Exact_circular_kernel_basic.cpp" )

else()
  
    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
  
endif()

