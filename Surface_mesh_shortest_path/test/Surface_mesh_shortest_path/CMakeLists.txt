# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Surface_mesh_shortest_path_ )

cmake_minimum_required(VERSION 2.8.11)


find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

    include_directories (BEFORE "include")

  create_single_source_cgal_program( "Surface_mesh_shortest_path_test_1.cpp" )
  create_single_source_cgal_program( "Surface_mesh_shortest_path_test_2.cpp" )
  create_single_source_cgal_program( "Surface_mesh_shortest_path_test_3.cpp" )
  create_single_source_cgal_program( "Surface_mesh_shortest_path_test_4.cpp" )
  create_single_source_cgal_program( "Surface_mesh_shortest_path_traits_test.cpp" )


  # Link with Boost.ProgramOptions (optional)
  find_package(Boost QUIET COMPONENTS program_options)
  if(Boost_PROGRAM_OPTIONS_FOUND)
    if( CGAL_AUTO_LINK_ENABLED )
      message( STATUS "Boost.ProgramOptions library: found" )
    else()
      message( STATUS "Boost.ProgramOptions library: ${Boost_PROGRAM_OPTIONS_LIBRARY}" )
    endif()
    add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
    list(APPEND CGAL_3RD_PARTY_LIBRARIES ${Boost_PROGRAM_OPTIONS_LIBRARY})
    create_single_source_cgal_program( "TestMesh.cpp" )
  else()
    message(STATUS "NOTICE: Example TestMesh.cpp requires boost program_option and will not be compiled.")
  endif()


else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

