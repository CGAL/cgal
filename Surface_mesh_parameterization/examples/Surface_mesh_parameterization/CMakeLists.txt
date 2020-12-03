# This is the CMake script for compiling this folder.

cmake_minimum_required(VERSION 3.1...3.15)
project(Surface_mesh_parameterization_Examples)

# Find CGAL
find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
include(CGAL_Eigen_support)
if(TARGET CGAL::Eigen_support)
  # Executables that require Eigen 3.1

  # ------------------------------------------------------------------
  # Detect SuiteSparse libraries:
  # If not found automatically, set SuiteSparse_DIR in CMake to the
  # directory where SuiteSparse was built.
  # ------------------------------------------------------------------

  set(SuiteSparse_USE_LAPACK_BLAS ON)
  find_package(SuiteSparse QUIET NO_MODULE
  )# 1st: Try to locate the *config.cmake file.
  if(NOT SuiteSparse_FOUND)
    set(SuiteSparse_VERBOSE ON)
    find_package(SuiteSparse QUIET) # 2nd: Use FindSuiteSparse.cmake module
    if(SuiteSparse_FOUND)
      include_directories(${SuiteSparse_INCLUDE_DIRS})
    endif(SuiteSparse_FOUND)
  else(NOT SuiteSparse_FOUND)
    include(${USE_SuiteSparse})
  endif(NOT SuiteSparse_FOUND)

  if(SuiteSparse_FOUND)
    if(SuiteSparse_UMFPACK_FOUND)
      message(STATUS "SuiteSparse_LIBS: ${SuiteSparse_LIBRARIES}")
      message(STATUS "Orbifold Tutte Embeddings will use UmfPackLU")
      add_definitions(-DEIGEN_DONT_ALIGN_STATICALLY)
      add_definitions(-DCGAL_SMP_USE_SPARSESUITE_SOLVERS)
    else()
      message(
        STATUS
          "NOTICE: The example `orbifold.cpp` will be compiled without the Sparsesuite library and UmfPack. Try setting SuiteSparse_UMF_INCLUDE_DIR and at least one of SuiteSparse_UMFPACK_LIBRARY_RELEASE and SuiteSparse_UMFPACK_LIBRARY_DEBUG to you UMFPACK installation."
      )
    endif()
  else(SuiteSparse_FOUND)
    message(
      STATUS
        "NOTICE: The example `orbifold.cpp` will be compiled without the Sparsesuite library."
    )
  endif(SuiteSparse_FOUND)

  # ------------------------------------------------------------------
  #   End of SuiteSparse detection
  # ------------------------------------------------------------------

  create_single_source_cgal_program("discrete_authalic.cpp")
  target_link_libraries(discrete_authalic PUBLIC CGAL::Eigen_support)
  create_single_source_cgal_program("lscm.cpp")
  target_link_libraries(lscm PUBLIC CGAL::Eigen_support)
  create_single_source_cgal_program("orbifold.cpp")
  target_link_libraries(orbifold PUBLIC CGAL::Eigen_support)
  create_single_source_cgal_program("seam_Polyhedron_3.cpp")
  target_link_libraries(seam_Polyhedron_3 PUBLIC CGAL::Eigen_support)
  create_single_source_cgal_program("simple_parameterization.cpp")
  target_link_libraries(simple_parameterization PUBLIC CGAL::Eigen_support)
  create_single_source_cgal_program("square_border_parameterizer.cpp")
  target_link_libraries(square_border_parameterizer PUBLIC CGAL::Eigen_support)
  create_single_source_cgal_program( "iterative_authalic_parameterizer.cpp" )
  target_link_libraries(iterative_authalic_parameterizer PUBLIC CGAL::Eigen_support)
  if(SuiteSparse_FOUND)
    target_link_libraries(orbifold PRIVATE ${SuiteSparse_LIBRARIES})
  endif()

else()
  message(
    STATUS
      "NOTICE: The examples require Eigen 3.1 (or greater) and will not be compiled."
  )
endif()
