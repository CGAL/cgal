if(EIGEN3_FOUND AND NOT TARGET CGAL::Eigen_support)
  if(NOT TARGET Threads::Threads)
    find_package(Threads REQUIRED)
  endif()
  add_library(CGAL::Eigen_support INTERFACE IMPORTED)
  set_target_properties(CGAL::Eigen_support PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "CGAL_EIGEN3_ENABLED"
    INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
endif()
