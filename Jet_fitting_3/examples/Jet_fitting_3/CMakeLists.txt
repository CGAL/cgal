# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Jet_fitting_3_example )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  # Link with BLAS and LAPACK (required)
  find_package(LAPACK)
  if(LAPACK_FOUND)

    include( ${LAPACK_USE_FILE} )

    # Link with Boost.ProgramOptions (optional)
    find_package(Boost QUIET COMPONENTS program_options)
    if(Boost_PROGRAM_OPTIONS_FOUND)
      if( CGAL_AUTO_LINK_ENABLED )
        message( STATUS "Boost.ProgramOptions library: found" )
      else()
        message( STATUS "Boost.ProgramOptions library: ${Boost_PROGRAM_OPTIONS_LIBRARY}" )
      endif()
      add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
      if ( NOT CGAL_AUTO_LINK_ENABLED )
        link_libraries( ${Boost_PROGRAM_OPTIONS_LIBRARY} )
      endif()
    endif()

    create_single_source_cgal_program( "Mesh_estimation.cpp" )
    create_single_source_cgal_program( "Single_estimation.cpp" )

  else(LAPACK_FOUND)

    message(STATUS "NOTICE: This program requires LAPACK, and will not be compiled.")

  endif(LAPACK_FOUND)

else()

    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

