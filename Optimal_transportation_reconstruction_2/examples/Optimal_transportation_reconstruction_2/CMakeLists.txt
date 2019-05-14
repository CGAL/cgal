project( Optimal_transportation_reconstruction_2_Examples )

cmake_minimum_required(VERSION 3.1)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  create_single_source_cgal_program( "otr2_simplest_example.cpp" )
  create_single_source_cgal_program( "otr2_simplest_example_with_tolerance.cpp" )
  create_single_source_cgal_program( "otr2_list_output_example.cpp" )
  create_single_source_cgal_program( "otr2_indexed_output_example.cpp" )
  create_single_source_cgal_program( "otr2_mass_example.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

