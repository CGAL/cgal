# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.15)
project( Surface_mesh_simplification_Tests )

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  find_package(Eigen3 3.2.0) #(requires 3.2.0 or greater)
  include(CGAL_Eigen_support)


include_directories( BEFORE  "c:/3rdPartyLibs/fast-envelope-VC/include" )
link_directories ( "c:/3rdPartyLibs/fast-envelope-VC/lib"  "c:/3rdPartyLibs/fast-envelope-VC"  "c:/3rdPartyLibs/fast-envelope-VC/tbb"   )

# create a target per cppfile
file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})
  create_single_source_cgal_program( "${cppfile}" )
endforeach()

  target_link_libraries( test_edge_collapse_FastEnvelope PUBLIC CGAL::Eigen_support)
  target_link_libraries( test_edge_collapse_FastEnvelope PUBLIC ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} FastEnvelope IndirectPredicates geogram)
else()

  message(STATUS "This program requires the CGAL library, and will not be compiled.")
  return()
endif()
