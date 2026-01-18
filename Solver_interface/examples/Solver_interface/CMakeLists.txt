# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Solver_interface_Examples)

find_package(CGAL REQUIRED)

# Use Eigen
find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("singular_value_decomposition.cpp")
  target_link_libraries(singular_value_decomposition PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("sparse_solvers.cpp")
  target_link_libraries(sparse_solvers PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("diagonalize_matrix.cpp")
  target_link_libraries(diagonalize_matrix PRIVATE CGAL::Eigen3_support)
else()
  message(STATUS "NOTICE: Eigen3 was not found. Some examples won't be available.")
endif()

find_package(OSQP QUIET)
include(CGAL_OSQP_support)
if(TARGET CGAL::OSQP_support)
  message(STATUS "Found OSQP")

  create_single_source_cgal_program("osqp_quadratic_program.cpp")
  target_link_libraries(osqp_quadratic_program PRIVATE CGAL::OSQP_support)
  message(STATUS "OSQP found and used")
else()
  message(STATUS "NOTICE: OSQP was not found. OSQP examples won't be available.")
endif()

find_package(SCIP QUIET)
include(CGAL_SCIP_support)
if(TARGET CGAL::SCIP_support)
  message(STATUS "Found SCIP")

  create_single_source_cgal_program("mixed_integer_program.cpp")
  target_link_libraries(mixed_integer_program PRIVATE CGAL::SCIP_support)
  message(STATUS "SCIP found and used")
else()
  find_package(GLPK QUIET)
  include(CGAL_GLPK_support)
  if(TARGET CGAL::GLPK_support)
    message(STATUS "Found GLPK")

    create_single_source_cgal_program("mixed_integer_program.cpp")
    target_link_libraries(mixed_integer_program PRIVATE CGAL::GLPK_support)
    message(STATUS "GLPK found and used")
  else()
    message(STATUS "NOTICE: This project requires either SCIP or GLPK, and will not be compiled. "
                   "Please provide either 'SCIP_DIR' or 'GLPK_INCLUDE_DIR' and 'GLPK_LIBRARIES'")
  endif()
endif()
