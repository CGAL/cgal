# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Polygon_mesh_processing_Tests )


cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

if ( COMMAND cmake_policy )

  cmake_policy( SET CMP0003 NEW )  

endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory

# include for local package
include_directories( BEFORE ../../include )

find_package(Eigen3 3.2.0) #(requires 3.2.0 or greater)

find_package( TBB )
if( TBB_FOUND )
  include( ${TBB_USE_FILE} )
  list( APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES} )
else()
  message( STATUS "NOTICE: Intel TBB was not found. test_pmp_distance will use sequential code." )
endif()

include( CGAL_CreateSingleSourceCGALProgram )
if (EIGEN3_FOUND)
  include( ${EIGEN3_USE_FILE} )

  # Creating entries for all .cpp/.C files with "main" routine
  # ##########################################################
  create_single_source_cgal_program("fairing_test.cpp")
  create_single_source_cgal_program("triangulate_hole_Polyhedron_3_no_delaunay_test.cpp" )
  create_single_source_cgal_program("triangulate_hole_Polyhedron_3_test.cpp")
endif(EIGEN3_FOUND)

  create_single_source_cgal_program("connected_component_polyhedron.cpp")
  create_single_source_cgal_program("connected_component_surface_mesh.cpp")
  create_single_source_cgal_program("pmp_compute_normals_test.cpp")
  create_single_source_cgal_program("orient_polygon_mesh_test.cpp" )
  create_single_source_cgal_program("orient_polygon_soup_test.cpp")
  create_single_source_cgal_program("point_inside_polyhedron_boundary_test.cpp")
  create_single_source_cgal_program("point_inside_surface_mesh_test.cpp" )
  create_single_source_cgal_program("polygon_mesh_slicer_test.cpp" )
  create_single_source_cgal_program("self_intersection_polyhedron_test.cpp" )
  create_single_source_cgal_program("self_intersection_surface_mesh_test.cpp" )
  create_single_source_cgal_program("test_is_polygon_soup_a_polygon_mesh.cpp")
  create_single_source_cgal_program("test_stitching.cpp")
  create_single_source_cgal_program("test_pmp_bgl_named_params.cpp")
  create_single_source_cgal_program("remeshing_test.cpp" )
  create_single_source_cgal_program("measures_test.cpp")
  create_single_source_cgal_program("triangulate_faces_test.cpp")
  create_single_source_cgal_program("test_pmp_remove_border_edge.cpp")
  create_single_source_cgal_program("test_pmp_distance.cpp")
  create_single_source_cgal_program("test_corefinement_and_constraints.cpp")
  create_single_source_cgal_program("test_corefinement_bool_op.cpp")
  create_single_source_cgal_program("test_corefine.cpp")
  create_single_source_cgal_program("test_coref_epic_points_identity.cpp")
  create_single_source_cgal_program("test_does_bound_a_volume.cpp")
  create_single_source_cgal_program("test_pmp_clip.cpp")
  create_single_source_cgal_program("triangulate_hole_polyline_test.cpp")
  create_single_source_cgal_program("test_comp_remeshing.cpp")
