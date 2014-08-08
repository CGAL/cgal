# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Voronoi_covariance_3_ )

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

  create_single_source_cgal_program( "convex_hull.cpp" )
  create_single_source_cgal_program( "convex_hull_translated.cpp" )
  create_single_source_cgal_program( "halfspaces_intersection.cpp" )
  create_single_source_cgal_program( "halfspaces_intersection_with_construction.cpp" )
  create_single_source_cgal_program( "lloyd_algorithm.cpp" )

  find_package(Boost COMPONENTS program_options filesystem graph)
  if (${Boost_PROGRAM_OPTIONS_FOUND} AND ${Boost_FILESYSTEM_FOUND})
    include_directories(${Boost_PROGRAM_OPTIONS_INCLUDE_DIR} ${Boost_FILESYSTEM_INCLUDE_DIR})
    create_single_source_cgal_program( "vcm.cpp" )
    target_link_libraries(vcm ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
  else()
    message(STATUS "NOTICE: An example needs Boost program options and filesystem and will not be compiled.")
  endif()

  find_package(Eigen3 3.1.0)
  if (EIGEN3_FOUND AND ${Boost_GRAPH_FOUND})
    include( ${EIGEN3_USE_FILE} )
    include_directories(${Boost_GRAPH_INCLUDE_DIR})
    create_single_source_cgal_program( "edges_estimation.cpp" )
    target_link_libraries(edges_estimation ${Boost_GRAPH_LIBRARY})
  else()
      message(STATUS "NOTICE: An example needs Eigen 3.1 (or greater) and BGL and will not be compiled.")
  endif()

else()
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()
