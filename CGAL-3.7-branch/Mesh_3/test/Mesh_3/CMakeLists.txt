# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Mesh_3_test )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

include_directories(../../include)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
find_package(CGAL COMPONENTS ImageIO)

if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )
  find_package(Boost)

  if ( Boost_FOUND AND Boost_VERSION GREATER 103400 )
    include( CGAL_CreateSingleSourceCGALProgram )

    # Test
    create_single_source_cgal_program("test_meshing_polyhedron.cpp")
    if( CGAL_ImageIO_USE_ZLIB )
      create_single_source_cgal_program("test_meshing_3D_image.cpp")
    else()
      message( STATUS "NOTICE: The example mesh_3D_image.cpp needs CGAL_ImageIO to be configured with ZLIB support, and will not be compiled." )
    endif()
    create_single_source_cgal_program("test_meshing_implicit_function.cpp")
    create_single_source_cgal_program("test_c3t3.cpp")
    create_single_source_cgal_program("test_criteria.cpp")
    create_single_source_cgal_program("test_robust_weighted_circumcenter.cpp")

  else()
    message(STATUS "NOTICE: This program requires Boost >= 1.34.1, and will not be compiled.")
  endif()

else()
  message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()

