# This is the CMake script for compiling this folder.


project( Surface_mesh_parameterization_Examples )

cmake_minimum_required(VERSION 2.8.11)


# Include this package's headers first
include_directories (BEFORE . include ../../include)

# Find CGAL
find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
  find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
 
  if (EIGEN3_FOUND)
    # Executables that require Eigen 3.1
    include( ${EIGEN3_USE_FILE} ) 

# ------------------------------------------------------------------
# Detect SuiteSparse libraries:
# If not found automatically, set SuiteSparse_DIR in CMake to the 
# directory where SuiteSparse was built.
# ------------------------------------------------------------------


set(SuiteSparse_USE_LAPACK_BLAS ON)
find_package(SuiteSparse QUIET NO_MODULE)  # 1st: Try to locate the *config.cmake file.
if(NOT SuiteSparse_FOUND)
        set(SuiteSparse_VERBOSE ON)
        find_package(SuiteSparse REQUIRED) # 2nd: Use FindSuiteSparse.cmake module
        include_directories(${SuiteSparse_INCLUDE_DIRS})
else()
        message(STATUS "Find SuiteSparse : include(${USE_SuiteSparse})")
        include(${USE_SuiteSparse})
endif()
MESSAGE(STATUS "SuiteSparse_LIBS: ${SuiteSparse_LIBRARIES}")
# ------------------------------------------------------------------
#   End of SuiteSparse detection
# ------------------------------------------------------------------





 
  create_single_source_cgal_program( "discrete_authalic.cpp" )
  create_single_source_cgal_program( "lscm.cpp" )
  create_single_source_cgal_program( "seam_Polyhedron_3.cpp" )
  create_single_source_cgal_program( "simple_parameterization.cpp" )
  create_single_source_cgal_program( "square_border_parameterizer.cpp" )
  create_single_source_cgal_program( "orbital.cpp" )
  
  target_link_libraries(orbital ${SuiteSparse_LIBRARIES} )
  
  else(EIGEN3_FOUND)
    message(STATUS "NOTICE: The examples require Eigen 3.1 (or greater) and will not be compiled.")  
  endif(EIGEN3_FOUND)

else()

  message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

