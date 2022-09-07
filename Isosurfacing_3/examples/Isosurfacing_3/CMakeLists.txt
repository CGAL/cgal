# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.14)
project( Isosurfacing_3_examples )

find_package(CGAL REQUIRED)

create_single_source_cgal_program( "marching_cubes_implicit_sphere.cpp" )
create_single_source_cgal_program( "marching_cubes_cartesian_grid_sphere.cpp" )
create_single_source_cgal_program( "marching_cubes_mesh_offset.cpp" )
create_single_source_cgal_program( "marching_cubes_inrimage.cpp" )
create_single_source_cgal_program( "dual_contouring_cartesian_grid.cpp" )
create_single_source_cgal_program( "dual_contouring_mesh_offset.cpp" )
create_single_source_cgal_program( "dual_contouring_octree.cpp" )
create_single_source_cgal_program( "all_cartesian_cube.cpp" )

find_package(TBB)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
  target_link_libraries(marching_cubes_implicit_sphere PUBLIC CGAL::TBB_support)
  target_link_libraries(marching_cubes_cartesian_grid_sphere PUBLIC CGAL::TBB_support)
  target_link_libraries(marching_cubes_mesh_offset PUBLIC CGAL::TBB_support)
  target_link_libraries(marching_cubes_inrimage PUBLIC CGAL::TBB_support)
  target_link_libraries(dual_contouring_cartesian_grid PUBLIC CGAL::TBB_support)
  target_link_libraries(dual_contouring_mesh_offset PUBLIC CGAL::TBB_support)
  target_link_libraries(dual_contouring_octree PUBLIC CGAL::TBB_support)
  target_link_libraries(all_cartesian_cube PUBLIC CGAL::TBB_support)
endif()
