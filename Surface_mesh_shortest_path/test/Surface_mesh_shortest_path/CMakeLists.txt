cmake_minimum_required(VERSION 3.1...3.13)

project( Surface_mesh_shortest_path_Tests )

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Core)

include(${CGAL_USE_FILE})

find_package(LEDA QUIET)

include_directories (BEFORE "include")

create_single_source_cgal_program( "Surface_mesh_shortest_path_test_1.cpp" )
create_single_source_cgal_program( "Surface_mesh_shortest_path_test_2.cpp" )
create_single_source_cgal_program( "Surface_mesh_shortest_path_test_3.cpp" )
create_single_source_cgal_program( "Surface_mesh_shortest_path_test_4.cpp" )
create_single_source_cgal_program( "Surface_mesh_shortest_path_traits_test.cpp" )

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
  if(CGAL_Core_FOUND OR LEDA_FOUND)
    create_single_source_cgal_program( "TestMesh.cpp" )
  else()
    message(STATUS "NOTICE: Example TestMesh.cpp requires the CGAL_Core library (or LEDA) and will not be compiled.")
  endif()
else()
  message(STATUS "NOTICE: Example TestMesh.cpp requires boost program_option and will not be compiled.")
endif()
