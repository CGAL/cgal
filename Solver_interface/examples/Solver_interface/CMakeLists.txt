# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.22)
project(Solver_interface_Examples)

find_package(CGAL REQUIRED)

# Use Eigen
find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
include(CGAL_Eigen3_support)

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("singular_value_decomposition.cpp")
  target_link_libraries(singular_value_decomposition PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("sparse_solvers.cpp")
  target_link_libraries(sparse_solvers PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("diagonalize_matrix.cpp")
  target_link_libraries(diagonalize_matrix PUBLIC CGAL::Eigen3_support)
endif()

find_package(OSQP QUIET)
include(CGAL_OSQP_support)

if(TARGET CGAL::OSQP_support)

  create_single_source_cgal_program("osqp_quadratic_program.cpp")
  target_link_libraries(osqp_quadratic_program PUBLIC CGAL::OSQP_support)
  message("OSQP found and used")

else()

  message(STATUS "NOTICE: OSQP was not found. OSQP examples won't be available.")

endif()

find_package(SCIP QUIET)
include(CGAL_SCIP_support)

if(TARGET CGAL::SCIP_support)

  create_single_source_cgal_program("mixed_integer_program.cpp")
  target_link_libraries(mixed_integer_program PUBLIC CGAL::SCIP_support)
  message("SCIP found and used")

else()

  find_package(GLPK QUIET)
  include(CGAL_GLPK_support)

  if(TARGET CGAL::GLPK_support)

    create_single_source_cgal_program("mixed_integer_program.cpp")
    target_link_libraries(mixed_integer_program PUBLIC CGAL::GLPK_support)
    message("GLPK found and used")

  else()

    message(
      STATUS
        "NOTICE : This project requires either SCIP or GLPK, and will not be compiled. "
        "Please provide either 'SCIP_DIR' or 'GLPK_INCLUDE_DIR' and 'GLPK_LIBRARIES'"
    )

  endif()

endif()
