# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Solver_interface_Examples )

cmake_minimum_required(VERSION 3.1)

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
  create_single_source_cgal_program( "mixed_integer_program.cpp" )

  find_package( SCIP )

  if (SCIP_FOUND)

    include_directories( BEFORE ${SCIP_INCLUDE_DIRS} )

    target_link_libraries( mixed_integer_program PRIVATE ${SCIP_LIBRARIES} )

    add_compile_definitions(CGAL_USE_SCIP)

    message("SCIP found and used")

  else()

      find_package( GLPK )

      if (GLPK_FOUND)

        include_directories( BEFORE ${GLPK_INCLUDE_DIR} )

        target_link_libraries( mixed_integer_program PRIVATE ${GLPK_LIBRARIES} )

        add_compile_definitions(CGAL_USE_GLPK)

        message("GLPK found and used")

      else()

        message(STATUS "This project requires either SCIP or GLPK, and will not be compiled. "
            "Please proivde either 'SCIP_DIR' or 'GLPK_INCLUDE_DIR' and 'GLPK_LIBRARIES'")

      endif()

  endif()

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

