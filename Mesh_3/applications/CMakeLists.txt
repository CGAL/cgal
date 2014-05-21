# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Mesh_3_applications )

CMAKE_MINIMUM_REQUIRED(VERSION 2.6.2)

include_directories(../include)


find_package(CGAL COMPONENTS ImageIO)

if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )

  find_package( Boost COMPONENTS program_options filesystem system)

  if ( Boost_FOUND AND Boost_VERSION GREATER 103400)
    if ( Boost_PROGRAM_OPTIONS_FOUND AND Boost_FILESYSTEM_FOUND )
      include_directories(${Boost_INCLUDE_DIRS})
      list(APPEND CGAL_3RD_PARTY_LIBRARIES "${Boost_LIBRARIES}")
    else()
      message(FATAL_ERROR "This program requires libs boost_program_options, boost_filesystem and boost_system")
    endif()
  else()
    message(FATAL_ERROR "This program requires Boost >= 1.34.1")
  endif()

  include( CGAL_CreateSingleSourceCGALProgram )

# Compilable examples
  create_single_source_cgal_program( "mesh_implicit_domains.cpp" "../examples/Mesh_3/implicit_functions.cpp" )
  if( CGAL_ImageIO_USE_ZLIB )
    create_single_source_cgal_program( "mesh_3D_image.cpp" )
  else()
    message( STATUS "NOTICE: The example mesh_3D_image.cpp needs CGAL_ImageIO to be configured with ZLIB support, and will not be compiled." )
  endif()
  create_single_source_cgal_program( "mesh_polyhedral_domain.cpp" )
  create_single_source_cgal_program( "output_distribution_to_stdout.cpp" )
  create_single_source_cgal_program( "mesher_tester.cpp" "../examples/Mesh_3/implicit_functions.cpp" )
  create_single_source_cgal_program( "mesher_tester_image.cpp" )
  create_single_source_cgal_program( "mesher_tester_polyhedron.cpp" )
else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()
