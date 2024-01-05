# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project( Isosurfacing_3_Examples )

find_package(CGAL REQUIRED)

create_single_source_cgal_program( "marching_cubes_implicit_sphere.cpp" )
create_single_source_cgal_program( "marching_cubes_seq_vs_parallel.cpp" )
create_single_source_cgal_program( "marching_cubes_Cartesian_grid_sphere.cpp" )
create_single_source_cgal_program( "marching_cubes_signed_mesh_offset.cpp" )
create_single_source_cgal_program( "marching_cubes_multiple_mesh_offsets.cpp" )
create_single_source_cgal_program( "marching_cubes_inrimage.cpp" )


find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program( "dual_contouring_Cartesian_grid.cpp" )
  target_link_libraries(dual_contouring_Cartesian_grid PRIVATE CGAL::Eigen3_support)

  create_single_source_cgal_program( "dual_contouring_mesh_offset.cpp" )
  target_link_libraries(dual_contouring_mesh_offset PRIVATE CGAL::Eigen3_support)

  create_single_source_cgal_program( "dual_contouring_octree.cpp" )
  target_link_libraries(dual_contouring_octree PRIVATE CGAL::Eigen3_support)

  create_single_source_cgal_program( "all_Cartesian_cube.cpp" )
  target_link_libraries(all_Cartesian_cube PRIVATE CGAL::Eigen3_support)

  create_single_source_cgal_program( "dual_contouring_implicit_iwp.cpp" )
  target_link_libraries(dual_contouring_implicit_iwp PRIVATE CGAL::Eigen3_support)

else()
  message(STATUS "NOTICE: Some examples use Eigen, and will not be compiled.")
endif()

find_package(TBB QUIET)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
  target_link_libraries(marching_cubes_implicit_sphere PRIVATE CGAL::TBB_support)
  target_link_libraries(marching_cubes_Cartesian_grid_sphere PRIVATE CGAL::TBB_support)
  target_link_libraries(marching_cubes_signed_mesh_offset PRIVATE CGAL::TBB_support)
  target_link_libraries(marching_cubes_inrimage PRIVATE CGAL::TBB_support)

  if(TARGET CGAL::Eigen3_support)
    target_link_libraries(dual_contouring_Cartesian_grid PRIVATE CGAL::TBB_support)
    target_link_libraries(dual_contouring_mesh_offset PRIVATE CGAL::TBB_support)
    target_link_libraries(dual_contouring_octree PRIVATE CGAL::TBB_support)
    target_link_libraries(all_Cartesian_cube PRIVATE CGAL::TBB_support)
    target_link_libraries(dual_contouring_implicit_iwp PRIVATE CGAL::TBB_support)
  endif()
endif()
