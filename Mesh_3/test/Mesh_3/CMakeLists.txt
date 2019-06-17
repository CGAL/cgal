# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Mesh_3_Tests )

cmake_minimum_required(VERSION 3.1)

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

find_package(CGAL QUIET COMPONENTS ImageIO)

if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )

  find_package( TBB QUIET )

  # Use Eigen
  find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
  if (NOT EIGEN3_FOUND)
    message(STATUS "This project requires the Eigen library, and will not be compiled.")
    return()
  else()
    include( ${EIGEN3_USE_FILE} )
  endif()
  
  create_single_source_cgal_program( "test_boost_has_xxx.cpp" )
  create_single_source_cgal_program( "test_c3t3.cpp" )
  create_single_source_cgal_program( "test_mesh_capsule_var_distance_bound.cpp" )
  create_single_source_cgal_program( "test_implicit_multi_domain_to_labeling_function_wrapper.cpp" )
  create_single_source_cgal_program( "test_c3t3_io.cpp" )
  create_single_source_cgal_program( "test_c3t3_with_features.cpp" )
  create_single_source_cgal_program( "test_criteria.cpp" )
  create_single_source_cgal_program( "test_domain_with_polyline_features.cpp" )
  create_single_source_cgal_program( "test_labeled_mesh_domain_3.cpp" )
  create_single_source_cgal_program( "test_mesh_criteria_creation.cpp" )
  create_single_source_cgal_program( "test_c3t3_into_facegraph.cpp" )
  create_single_source_cgal_program( "test_without_detect_features.cpp" )
  if(CGAL_ImageIO_USE_ZLIB)
    create_single_source_cgal_program( "test_meshing_3D_image.cpp" )
    create_single_source_cgal_program( "test_meshing_3D_image_deprecated.cpp" )
    create_single_source_cgal_program( "test_meshing_3D_gray_image.cpp" )
    create_single_source_cgal_program( "test_meshing_3D_gray_image_deprecated.cpp" )
  else()
    message(STATUS "test_meshing_3D_image requires the ZLIB library, and will not be compiled.")
  endif()
  create_single_source_cgal_program( "test_meshing_implicit_function.cpp" )
  create_single_source_cgal_program( "test_meshing_implicit_function_deprecated.cpp" )
  create_single_source_cgal_program( "test_meshing_polyhedral_complex.cpp" )
  create_single_source_cgal_program( "test_meshing_polyhedron.cpp" )
  create_single_source_cgal_program( "test_meshing_polylines_only.cpp" )
  create_single_source_cgal_program( "test_meshing_polyhedron_with_features.cpp" )
  create_single_source_cgal_program( "test_meshing_verbose.cpp" )
  create_single_source_cgal_program( "test_meshing_unit_tetrahedron.cpp" )
  create_single_source_cgal_program( "test_meshing_with_default_edge_size.cpp" )
  create_single_source_cgal_program( "test_meshing_determinism.cpp" )
  create_single_source_cgal_program( "test_c3t3_extract_subdomains_boundaries.cpp" )
  create_single_source_cgal_program( "test_mesh_3_issue_1554.cpp" )
  create_single_source_cgal_program( "test_mesh_polyhedral_domain_with_features_deprecated.cpp" )
  create_single_source_cgal_program( "test_meshing_with_one_step.cpp" )

  foreach(target
      test_meshing_verbose
      test_meshing_polyhedron_with_features
      test_meshing_utilities.h
      test_mesh_implicit_domains
      test_meshing_implicit_function
      test_meshing_3D_image
      test_meshing_3D_gray_image
      test_meshing_unit_tetrahedron
      test_backward_compatibility
      test_meshing_polyhedron
      test_meshing_polyhedral_complex
      )
    if(TBB_FOUND AND TARGET ${target})
      CGAL_target_use_TBB(${target})
    endif()
  endforeach()

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

