# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.15)
project( Solver_interface_Examples )


find_package(CGAL QUIET)

if ( CGAL_FOUND )

  # Use Eigen
  find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)

  if (EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )
    create_single_source_cgal_program( "singular_value_decomposition.cpp" )
    create_single_source_cgal_program( "sparse_solvers.cpp" )
  endif()

  create_single_source_cgal_program( "diagonalize_matrix.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

