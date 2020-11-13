cmake_minimum_required(VERSION 3.1...3.15)
project(Scale_space_reconstruction_3_Examples)

find_package(CGAL REQUIRED)

# Activate concurrency?
option(ACTIVATE_CONCURRENCY "Enable concurrency" ON)

if(ACTIVATE_CONCURRENCY)
  find_package(TBB)
  include(CGAL_TBB_support)
  if(NOT TARGET CGAL::TBB_support)
    message(
      STATUS
        "NOTICE: Intel TBB NOT found! The example is faster if TBB is linked.")
  endif()
endif()

find_package(Eigen3 3.1.0)
include(CGAL_Eigen_support)
if(TARGET CGAL::Eigen_support)
  create_single_source_cgal_program("scale_space.cpp")
  target_link_libraries(scale_space PUBLIC CGAL::Eigen_support)
  create_single_source_cgal_program("scale_space_incremental.cpp")
  target_link_libraries(scale_space_incremental PUBLIC CGAL::Eigen_support)
  create_single_source_cgal_program("scale_space_manifold.cpp")
  target_link_libraries(scale_space_manifold PUBLIC CGAL::Eigen_support)
  create_single_source_cgal_program("scale_space_advancing_front.cpp")
  target_link_libraries(scale_space_advancing_front PUBLIC CGAL::Eigen_support)
  if(ACTIVATE_CONCURRENCY AND TARGET CGAL::TBB_support)
    target_link_libraries(scale_space PUBLIC CGAL::TBB_support)
    target_link_libraries(scale_space_incremental PUBLIC CGAL::TBB_support)
    target_link_libraries(scale_space_manifold PUBLIC CGAL::TBB_support)
    target_link_libraries(scale_space_advancing_front PUBLIC CGAL::TBB_support)
  endif()
else()
  message(
    STATUS
      "NOTICE: The example needs Eigen 3.1 (or greater) and will not be compiled."
  )
endif()
