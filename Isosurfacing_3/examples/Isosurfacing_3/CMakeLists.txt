# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project( Isosurfacing_3_Examples )

find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)

find_package(TBB QUIET)
include(CGAL_TBB_support)

create_single_source_cgal_program("marching_cubes.cpp")

# undocumented
create_single_source_cgal_program("marching_cubes_strategies.cpp")

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("dual_contouring.cpp")
  create_single_source_cgal_program("contouring_discrete_data.cpp")
  create_single_source_cgal_program("contouring_image.cpp")
  create_single_source_cgal_program("contouring_implicit_data.cpp")
  create_single_source_cgal_program("contouring_mesh_offset.cpp")
  create_single_source_cgal_program("contouring_octree.cpp")

  # undocumented
  create_single_source_cgal_program("dual_contouring_strategies.cpp")
  create_single_source_cgal_program("dual_contouring_intersection_oracles.cpp")

  target_link_libraries(dual_contouring PRIVATE CGAL::Eigen3_support)
  target_link_libraries(contouring_discrete_data PRIVATE CGAL::Eigen3_support)
  target_link_libraries(contouring_image PRIVATE CGAL::Eigen3_support)
  target_link_libraries(contouring_implicit_data PRIVATE CGAL::Eigen3_support)
  target_link_libraries(contouring_mesh_offset PRIVATE CGAL::Eigen3_support)
  target_link_libraries(contouring_octree PRIVATE CGAL::Eigen3_support)

  target_link_libraries(dual_contouring_strategies PRIVATE CGAL::Eigen3_support)
  target_link_libraries(dual_contouring_intersection_oracles PRIVATE CGAL::Eigen3_support)

  if(TARGET CGAL::TBB_support)
    target_link_libraries(dual_contouring PRIVATE CGAL::TBB_support)
    target_link_libraries(contouring_discrete_data PRIVATE CGAL::TBB_support)
    target_link_libraries(contouring_image PRIVATE CGAL::TBB_support)
    target_link_libraries(contouring_implicit_data PRIVATE CGAL::TBB_support)
    target_link_libraries(contouring_mesh_offset PRIVATE CGAL::TBB_support)
    target_link_libraries(contouring_octree PRIVATE CGAL::TBB_support)

    target_link_libraries(dual_contouring_strategies PRIVATE CGAL::TBB_support)
    target_link_libraries(dual_contouring_intersection_oracles PRIVATE CGAL::TBB_support)
  endif()
else()
  message(STATUS "NOTICE: Some examples use Eigen, and will not be compiled.")
endif()

if(TARGET CGAL::TBB_support)
  target_link_libraries(marching_cubes PRIVATE CGAL::TBB_support)
  target_link_libraries(marching_cubes_strategies PRIVATE CGAL::TBB_support)
endif()
