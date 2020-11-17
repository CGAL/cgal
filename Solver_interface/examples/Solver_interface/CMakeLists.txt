# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Solver_interface_Examples)

find_package(CGAL REQUIRED)

# Use Eigen
find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
include(CGAL_Eigen_support)

if(TARGET CGAL::Eigen_support)
  create_single_source_cgal_program("singular_value_decomposition.cpp")
  target_link_libraries(singular_value_decomposition PUBLIC CGAL::Eigen_support)
  create_single_source_cgal_program("sparse_solvers.cpp")
  target_link_libraries(sparse_solvers PUBLIC CGAL::Eigen_support)
  create_single_source_cgal_program("diagonalize_matrix.cpp")
  target_link_libraries(diagonalize_matrix PUBLIC CGAL::Eigen_support)
endif()

create_single_source_cgal_program("mixed_integer_program.cpp")

find_package(SCIP QUIET)
include(CGAL_SCIP_support)

if(TARGET CGAL::SCIP_support)

  target_link_libraries(mixed_integer_program PUBLIC CGAL::SCIP_support)
  message("SCIP found and used")

else()

  find_package(GLPK QUIET)
  include(CGAL_GLPK_support)

  if(TARGET CGAL::GLPK_support)

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
