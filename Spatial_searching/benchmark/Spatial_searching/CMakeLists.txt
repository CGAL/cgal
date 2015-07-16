# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Spatial_searching_ )

cmake_minimum_required(VERSION 2.8.11)


find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
          
find_package(Eigen3 3.1.91) #(requires 3.2.0 or greater)
if (EIGEN3_FOUND)
  include( ${EIGEN3_USE_FILE} )
endif()

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include" "./include")

    include_directories (BEFORE "include")

  create_single_source_cgal_program( "Compare_ANN_STANN_CGAL.cpp" )
  create_single_source_cgal_program( "nanoflan.cpp" )
create_single_source_cgal_program( "binary.cpp" )
  create_single_source_cgal_program( "nearest_neighbor_searching_50.cpp" )
  create_single_source_cgal_program( "nearest_neighbor_searching_inplace_50.cpp" )
  create_single_source_cgal_program( "Nearest_neighbor_searching.cpp" )
  create_single_source_cgal_program( "Nearest_neighbor_searching_2D.cpp" )
  create_single_source_cgal_program( "Nearest_neighbor_searching_2D_user_defined.cpp" )
  create_single_source_cgal_program( "Split_data.cpp" )
create_single_source_cgal_program( "nn3cgal.cpp" )
create_single_source_cgal_program( "nn4cgal.cpp" )
create_single_source_cgal_program( "nn3nanoflan.cpp" )
create_single_source_cgal_program( "sizeof.cpp" )
create_single_source_cgal_program( "deque.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

