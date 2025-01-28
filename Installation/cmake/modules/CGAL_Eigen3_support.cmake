if((EIGEN3_FOUND OR Eigen3_FOUND) AND NOT TARGET CGAL::Eigen3_support)
  if(NOT TARGET Threads::Threads)
    find_package(Threads REQUIRED)
  endif()
  add_library(CGAL::Eigen3_support INTERFACE IMPORTED)
  set_target_properties(CGAL::Eigen3_support PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "CGAL_EIGEN3_ENABLED")
  if(TARGET Eigen3::Eigen)
    target_link_libraries(CGAL::Eigen3_support INTERFACE Eigen3::Eigen)
  else()
    set_target_properties(CGAL::Eigen3_support PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
  endif()
endif()
