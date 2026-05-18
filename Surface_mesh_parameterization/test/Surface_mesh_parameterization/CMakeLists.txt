# This is the CMake script for compiling this folder.

cmake_minimum_required(VERSION 3.12...3.31)
project(Surface_mesh_parameterization_Tests)

# Find CGAL
find_package(CGAL REQUIRED)

find_package(UMFPACK QUIET NO_MODULE)

if(TARGET SuiteSparse::UMFPACK)
  message(STATUS "Orbifold Tutte Embeddings will use UmfPackLU")
else()
  message(STATUS "NOTICE: Examples will be compiled without the SuiteSparse library and UmfPack. Try setting CMAKE_PREFIX_PATH to your SuiteSparse installation.")
endif()

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("extensive_parameterization_test.cpp")
  target_link_libraries(extensive_parameterization_test PRIVATE CGAL::Eigen3_support)

  if(TARGET SuiteSparse::UMFPACK)
    target_link_libraries(extensive_parameterization_test PRIVATE SuiteSparse::UMFPACK)
    target_compile_definitions(extensive_parameterization_test PRIVATE EIGEN_DONT_ALIGN_STATICALLY CGAL_SMP_USE_SUITESPARSE_SOLVERS)
  endif()
else()
  message("NOTICE: The tests require Eigen 3.1 (or greater), and will not be compiled.")
endif()
