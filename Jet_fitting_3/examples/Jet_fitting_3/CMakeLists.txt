# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Jet_fitting_3_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CreateSingleSourceCGALProgram )
  include( GeneratorSpecificSettings )

  # Link with BLAS and LAPACK (required)
  include( CGAL_UseLAPACK )

  # Link with Boost.ProgramOptions (optional)
  include( CGAL_UseBoostProgramOptions )

  create_single_source_cgal_program( "Mesh_estimation.cpp" )
  create_single_source_cgal_program( "Single_estimation.cpp" )

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

