# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Segment_Delaunay_graph_Linf_2_Examples )

cmake_minimum_required(VERSION 3.1)


find_package(CGAL QUIET)

if ( CGAL_FOUND )

  create_single_source_cgal_program( "sdg-count-sites-linf.cpp" )
  create_single_source_cgal_program( "sdg-fast-sp-linf.cpp" )
  create_single_source_cgal_program( "sdg-fast-sp-polygon-linf.cpp" )
  create_single_source_cgal_program( "sdg-filtered-traits-linf.cpp" )
  create_single_source_cgal_program( "sdg-info-set-linf.cpp" )
  create_single_source_cgal_program( "sdg-red-blue-info-linf.cpp" )
  create_single_source_cgal_program( "sdg-voronoi-edges-linf.cpp" )
  create_single_source_cgal_program( "sdg-voronoi-edges-exact-linf.cpp" )
  create_single_source_cgal_program( "print-sdg-linf.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

