# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.15)
project( Jet_fitting_3_Examples )



find_package(CGAL QUIET)

if ( CGAL_FOUND )

  # use Eigen
  find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
  if (EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )
    # Link with Boost.ProgramOptions (optional)
    find_package(Boost QUIET COMPONENTS program_options)
    if(Boost_PROGRAM_OPTIONS_FOUND)
      if(TARGET Boost::program_options)
        set(Boost_PROGRAM_OPTIONS_LIBRARY Boost::program_options)
      endif()
      if( CGAL_AUTO_LINK_ENABLED )
        message( STATUS "Boost.ProgramOptions library: found" )
      else()
        message( STATUS "Boost.ProgramOptions library: ${Boost_PROGRAM_OPTIONS_LIBRARY}" )
      endif()
      add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
      list(APPEND CGAL_3RD_PARTY_LIBRARIES ${Boost_PROGRAM_OPTIONS_LIBRARY})
   endif()

    create_single_source_cgal_program( "Mesh_estimation.cpp" )
    create_single_source_cgal_program( "Single_estimation.cpp" )

  else()
    message(STATUS "NOTICE: This program requires Eigen 3.1 (or greater) and will not be compiled.")
  endif()

else()

  message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

