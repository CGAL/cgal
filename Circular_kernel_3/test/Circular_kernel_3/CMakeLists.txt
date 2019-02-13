# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Circular_kernel_3_Tests )

CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

include(${CGAL_USE_FILE})

include_directories (BEFORE include)
include_directories (BEFORE ../../../Kernel_23/test/Kernel_23/include)

create_single_source_cgal_program( "test_Spherical_kernel.cpp" )
create_single_source_cgal_program( "test_Spherical_kernel_basics.cpp" )
create_single_source_cgal_program( "test_Exact_spherical_kernel.cpp" )
create_single_source_cgal_program( "test_Lazy_Spherical_kernel.cpp" )
create_single_source_cgal_program( "test_Lazy_spherical_kernel_basics.cpp" )
create_single_source_cgal_program( "test_Spherical_kernel_with_core.cpp" )

else()
  
    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
  
endif()

