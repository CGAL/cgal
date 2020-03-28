if (CGAL_target_use_Eigen_included)
  return()
endif()
set(CGAL_target_use_Eigen_included TRUE)

set( Eigen3_FIND_VERSION "3.1.0")
set(EIGEN3_USE_FILE "UseEigen3")

function(CGAL_target_use_Eigen target)
  target_include_directories(${target} PUBLIC ${EIGEN3_INCLUDE_DIR})
  target_compile_options( ${target} PUBLIC -DCGAL_EIGEN3_ENABLED)
endfunction()

