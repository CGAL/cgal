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
  include( CGAL_SetupBLAS )
  include( CGAL_SetupLAPACK )

  # Link with Boost.ProgramOptions (optional)
  find_package(Boost 1.33.1 QUIET COMPONENTS program_options)
  if(Boost_PROGRAM_OPTIONS_FOUND)
    add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
    set( CGAL_3RD_PARTY_DEFINITIONS  ${CGAL_3RD_PARTY_DEFINITIONS}  "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )

    if( AUTO_LINK_ENABLED )
      message( STATUS "Boost.ProgramOptions library: found" )
    else()
      message( STATUS "Boost.ProgramOptions library: ${Boost_PROGRAM_OPTIONS_LIBRARY}" )
      set( CGAL_3RD_PARTY_LIBRARIES    ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_PROGRAM_OPTIONS_LIBRARY} )
    endif()
  endif()

  create_single_source_cgal_program( "Mesh_estimation.cpp" )
  create_single_source_cgal_program( "Single_estimation.cpp" )

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

