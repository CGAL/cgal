# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Surface_mesh_skeletonization_ )

cmake_minimum_required(VERSION 2.8.11)


find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  find_package(Eigen3 3.2.0) #(requires 3.2.0 or greater)

  if(EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )

    include( CGAL_CreateSingleSourceCGALProgram )

    include_directories (BEFORE "../../include")

    create_single_source_cgal_program( "simple_mcfskel_example.cpp" )
    create_single_source_cgal_program( "simple_mcfskel_sm_example.cpp" )
    create_single_source_cgal_program( "MCF_Skeleton_example.cpp" )
    create_single_source_cgal_program( "MCF_Skeleton_sm_example.cpp" )
    create_single_source_cgal_program( "segmentation_example.cpp" )
  else()
    message(STATUS "These programs require the Eigen library (3.2 or greater), and will not be compiled.")
  endif()

else()

  message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

