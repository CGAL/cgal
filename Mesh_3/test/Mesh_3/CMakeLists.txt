# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.20)
project( Mesh_3_Tests )


if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

find_package(CGAL QUIET COMPONENTS ImageIO)

if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )

  find_package( TBB QUIET )
  include(CGAL_TBB_support)

  find_package( ITT QUIET )

  # Use Eigen
  find_package(Eigen3 3.1.0 REQUIRED) #(3.1.0 or greater)
  include(CGAL_Eigen3_support)
  if (NOT TARGET CGAL::Eigen3_support)
    message(STATUS "This project requires the Eigen library, and will not be compiled.")
    return()
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
  create_single_source_cgal_program( "test_mesh_cell_base_3.cpp")

  foreach(target
      test_boost_has_xxx
      test_c3t3
      test_mesh_capsule_var_distance_bound
      test_implicit_multi_domain_to_labeling_function_wrapper
      test_c3t3_io
      test_c3t3_with_features
      test_criteria
      test_domain_with_polyline_features
      test_labeled_mesh_domain_3
      test_mesh_criteria_creation
      test_c3t3_into_facegraph
      test_without_detect_features
      test_meshing_3D_image
      test_meshing_3D_image_deprecated
      test_meshing_3D_gray_image
      test_meshing_3D_gray_image_deprecated
      test_meshing_implicit_function
      test_meshing_implicit_function_deprecated
      test_meshing_polyhedral_complex
      test_meshing_polyhedron
      test_meshing_polylines_only
      test_meshing_polyhedron_with_features
      test_meshing_verbose
      test_meshing_unit_tetrahedron
      test_meshing_with_default_edge_size
      test_meshing_determinism
      test_c3t3_extract_subdomains_boundaries
      test_mesh_3_issue_1554
      test_mesh_polyhedral_domain_with_features_deprecated
      test_mesh_cell_base_3
      test_meshing_with_one_step.cpp)
    if(TARGET ${target})
      target_link_libraries(${target} PUBLIC CGAL::Eigen3_support)
    endif()
  endforeach()

  if(TARGET CGAL::TBB_support)
    foreach(target
        test_meshing_verbose
        test_meshing_polyhedron_with_features
        test_meshing_utilities.h
        test_meshing_implicit_function
        test_meshing_3D_image
        test_meshing_3D_gray_image
        test_meshing_unit_tetrahedron
        test_meshing_polyhedron
        test_meshing_polyhedral_complex
        test_mesh_capsule_var_distance_bound
        test_mesh_3_issue_1554
        test_mesh_polyhedral_domain_with_features_deprecated
        test_mesh_cell_base_3
        )
      if(TARGET ${target})
        target_link_libraries(${target} PUBLIC CGAL::TBB_support)
      endif()
    endforeach()

    if(BUILD_TESTING)
      set_property(TEST
        execution___of__test_meshing_verbose
        execution___of__test_meshing_polyhedron_with_features
        execution___of__test_meshing_implicit_function
        execution___of__test_meshing_3D_image
        execution___of__test_meshing_3D_gray_image
        execution___of__test_meshing_unit_tetrahedron
        execution___of__test_meshing_polyhedron
        execution___of__test_meshing_polyhedral_complex
        execution___of__test_mesh_capsule_var_distance_bound
        execution___of__test_mesh_3_issue_1554
        execution___of__test_mesh_polyhedral_domain_with_features_deprecated
        execution___of__test_mesh_cell_base_3
      PROPERTY RUN_SERIAL 1)
    endif()
  endif()
  if(TARGET ITT::ITT)
    target_link_libraries(test_meshing_polyhedron_with_features PRIVATE ITT::ITT)
    target_compile_definitions(test_meshing_polyhedron_with_features PRIVATE CGAL_MESH_3_USE_INTEL_ITT)
    target_link_libraries(test_meshing_verbose PRIVATE ITT::ITT)
    target_compile_definitions(test_meshing_verbose PRIVATE CGAL_MESH_3_USE_INTEL_ITT)
  endif()

  if(BUILD_TESTING)
    set_tests_properties(
      execution___of__test_meshing_polyhedron_with_features
      execution___of__test_meshing_verbose
      PROPERTIES RESOURCE_LOCK Mesh_3_Tests_IO)
  endif()
else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
