# This is the CMake script for compiling this folder.

cmake_minimum_required(VERSION 3.12...3.31)
project(Surface_mesh_parameterization_Examples)

find_package(CGAL REQUIRED)

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(TARGET CGAL::Eigen3_support)
  # Executables that require Eigen 3.1

  find_package(UMFPACK QUIET NO_MODULE)

  if(TARGET SuiteSparse::UMFPACK)
    message(STATUS "Orbifold Tutte Embeddings will use UmfPackLU")
  else()
    message(STATUS "NOTICE: Examples will be compiled without the SuiteSparse library and UmfPack. Try setting CMAKE_PREFIX_PATH to your SuiteSparse installation.")
  endif()

  create_single_source_cgal_program("ARAP_parameterization.cpp")
  target_link_libraries(ARAP_parameterization PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("discrete_authalic.cpp")
  target_link_libraries(discrete_authalic PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("lscm.cpp")
  target_link_libraries(lscm PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("orbifold.cpp")
  target_link_libraries(orbifold PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("seam_Polyhedron_3.cpp")
  target_link_libraries(seam_Polyhedron_3 PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("simple_parameterization.cpp")
  target_link_libraries(simple_parameterization PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program("square_border_parameterizer.cpp")
  target_link_libraries(square_border_parameterizer PRIVATE CGAL::Eigen3_support)
  create_single_source_cgal_program( "iterative_authalic_parameterizer.cpp" )
  target_link_libraries(iterative_authalic_parameterizer PRIVATE CGAL::Eigen3_support)
  if(TARGET SuiteSparse::UMFPACK)
    target_link_libraries(ARAP_parameterization PRIVATE SuiteSparse::UMFPACK)
    target_link_libraries(orbifold PRIVATE PRIVATE SuiteSparse::UMFPACK)
    target_compile_definitions(ARAP_parameterization PRIVATE EIGEN_DONT_ALIGN_STATICALLY CGAL_SMP_USE_SUITESPARSE_SOLVERS)
    target_compile_definitions(orbifold PRIVATE EIGEN_DONT_ALIGN_STATICALLY CGAL_SMP_USE_SUITESPARSE_SOLVERS)
  endif()

else()
  message("NOTICE: The examples require Eigen 3.1 (or greater), and will not be compiled.")
endif()
