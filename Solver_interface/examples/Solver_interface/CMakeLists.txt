# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.15)
project( Solver_interface_Examples )


find_package(CGAL QUIET)

if ( CGAL_FOUND )

  # Use Eigen
  find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)

  if (EIGEN3_FOUND)
    create_single_source_cgal_program( "singular_value_decomposition.cpp" )
    CGAL_target_use_Eigen(singular_value_decomposition)
    create_single_source_cgal_program( "sparse_solvers.cpp" )
    CGAL_target_use_Eigen(sparse_solvers)
    create_single_source_cgal_program( "diagonalize_matrix.cpp" )
    CGAL_target_use_Eigen(diagonalize_matrix)
  endif()

  create_single_source_cgal_program( "mixed_integer_program.cpp" )

  find_package( SCIP QUIET)

  if (SCIP_FOUND)

    CGAL_target_use_SCIP(mixed_integer_program)
    message("SCIP found and used")

  else()

      find_package( GLPK QUIET)

      if (GLPK_FOUND)
        
        CGAL_target_use_GLPK(mixed_integer_program)
        message("GLPK found and used")

      else()

        message(STATUS "NOTICE : This project requires either SCIP or GLPK, and will not be compiled. "
            "Please provide either 'SCIP_DIR' or 'GLPK_INCLUDE_DIR' and 'GLPK_LIBRARIES'")

      endif()

  endif()

else()
  
    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
  
endif()

