# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Arrangement_on_surface_2_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  # This package has problems with g++ strict-aliasing optimizations
  # (implied by -O2, -O3 and -Os).
  if (CMAKE_COMPILER_IS_GNUCXX)
    SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-strict-aliasing")
    SET (CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   -fno-strict-aliasing")
  endif()

  create_single_source_cgal_program( "aggregated_insertion.cpp" )
  create_single_source_cgal_program( "batched_point_location.cpp" )
  create_single_source_cgal_program( "Bezier_curves.cpp" )
  create_single_source_cgal_program( "bgl_dual_adapter.cpp" )
  create_single_source_cgal_program( "bgl_primal_adapter.cpp" )
  create_single_source_cgal_program( "bounded_planar_vertical_decomposition.cpp" )
  create_single_source_cgal_program( "circles.cpp" )
  create_single_source_cgal_program( "circular_arcs.cpp" )
  create_single_source_cgal_program( "circular_line_arcs.cpp" )
  create_single_source_cgal_program( "conic_multiplicities.cpp" )
  create_single_source_cgal_program( "conics.cpp" )
  create_single_source_cgal_program( "consolidated_curve_data.cpp" )
  create_single_source_cgal_program( "curve_history.cpp" )
  create_single_source_cgal_program( "dcel_extension.cpp" )
  create_single_source_cgal_program( "dcel_extension_io.cpp" )
  create_single_source_cgal_program( "dual_lines.cpp" )
  create_single_source_cgal_program( "dual_with_data.cpp" )
  create_single_source_cgal_program( "edge_insertion.cpp" )
  create_single_source_cgal_program( "edge_manipulation.cpp" )
  create_single_source_cgal_program( "edge_manipulation_curve_history.cpp" )
  create_single_source_cgal_program( "face_extension.cpp" )
  create_single_source_cgal_program( "face_extension_overlay.cpp" )
  create_single_source_cgal_program( "generic_curve_data.cpp" )
  create_single_source_cgal_program( "global_insertion.cpp" )
  create_single_source_cgal_program( "global_removal.cpp" )
  create_single_source_cgal_program( "incremental_insertion.cpp" )
  create_single_source_cgal_program( "io.cpp" )
  create_single_source_cgal_program( "io_curve_history.cpp" )
  create_single_source_cgal_program( "io_unbounded.cpp" )
  create_single_source_cgal_program( "isolated_vertices.cpp" )
  create_single_source_cgal_program( "observer.cpp" )
  create_single_source_cgal_program( "overlay.cpp" )
  create_single_source_cgal_program( "overlay_unbounded.cpp" )
  create_single_source_cgal_program( "point_location.cpp" )
  create_single_source_cgal_program( "polylines.cpp" )
  create_single_source_cgal_program( "predefined_kernel.cpp" )
  create_single_source_cgal_program( "predefined_kernel_non_intersecting.cpp" )
  create_single_source_cgal_program( "rational_functions.cpp" )
  create_single_source_cgal_program( "special_edge_insertion.cpp" )
  create_single_source_cgal_program( "spherical_insert.cpp" )
  create_single_source_cgal_program( "sweep_line.cpp" )
  create_single_source_cgal_program( "tracing_counting.cpp" )
  create_single_source_cgal_program( "unbounded_non_intersecting.cpp" )
  create_single_source_cgal_program( "unbounded_rational_functions.cpp" )
  create_single_source_cgal_program( "unb_planar_vertical_decomposition.cpp" )
  create_single_source_cgal_program( "vertical_ray_shooting.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

