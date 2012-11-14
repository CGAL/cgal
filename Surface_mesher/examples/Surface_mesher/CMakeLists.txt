# This is the CMake script for compiling a CGAL application.


project( Surface_mesher_example )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

include_directories(../../include)

find_package(CGAL QUIET COMPONENTS ImageIO)

if ( CGAL_FOUND AND CGAL_ImageIO_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  create_single_source_cgal_program( "mesh_a_3d_gray_image.cpp" )
  create_single_source_cgal_program( "mesh_an_implicit_function.cpp" )

else()
  if(RUNNING_CGAL_AUTO_TEST)
    # Just to avoid a warning from CMake if that variable is set on the command line...
  endif()
  
  message(STATUS "NOTICE: This program requires the CGAL and CGAL ImageIO libraries, and will not be compiled.")
  
endif()

