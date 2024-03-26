# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project( Mesh_3_Tests )

find_package(CGAL REQUIRED COMPONENTS ImageIO)

# Use Eigen
find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)
if (NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: This project requires the Eigen library, and will not be compiled.")
  return()
endif()

find_package(TBB QUIET)
include(CGAL_TBB_support)

find_package(ITT QUIET)

create_single_source_cgal_program( "test_boost_has_xxx.cpp" )
create_single_source_cgal_program( "test_mesh_capsule_var_distance_bound.cpp" )
create_single_source_cgal_program( "test_implicit_multi_domain_to_labeling_function_wrapper.cpp" )
create_single_source_cgal_program( "test_criteria.cpp" )
create_single_source_cgal_program( "test_domain_with_polyline_features.cpp" )
create_single_source_cgal_program( "test_labeled_mesh_domain_3.cpp" )
create_single_source_cgal_program( "test_mesh_criteria_creation.cpp" )
create_single_source_cgal_program( "test_without_detect_features.cpp" )

if(CGAL_ImageIO_USE_ZLIB)
  create_single_source_cgal_program( "test_meshing_3D_image.cpp" )
  create_single_source_cgal_program( "test_meshing_3D_image_with_features.cpp" )
  create_single_source_cgal_program( "test_meshing_3D_gray_image.cpp" )
  create_single_source_cgal_program( "test_min_size_criteria.cpp")
else()
  message(STATUS "NOTICE: The test 'test_meshing_3D_image' requires the ZLIB library, and will not be compiled.")
endif()

create_single_source_cgal_program( "test_meshing_implicit_function.cpp" )
create_single_source_cgal_program( "test_meshing_polyhedral_complex.cpp" )
create_single_source_cgal_program( "test_meshing_polyhedral_complex_with_manifold_and_min_size.cpp")
create_single_source_cgal_program( "test_meshing_polyhedron.cpp" )
create_single_source_cgal_program( "test_meshing_polylines_only.cpp" )
create_single_source_cgal_program( "test_meshing_polyhedron_with_features.cpp" )
create_single_source_cgal_program( "test_meshing_verbose.cpp" )
create_single_source_cgal_program( "test_meshing_unit_tetrahedron.cpp" )
create_single_source_cgal_program( "test_meshing_with_default_edge_size.cpp" )
create_single_source_cgal_program( "test_meshing_determinism.cpp" )
create_single_source_cgal_program( "test_meshing_without_features_determinism.cpp" )
create_single_source_cgal_program( "test_mesh_3_issue_1554.cpp" )
create_single_source_cgal_program( "test_mesh_polyhedral_domain_with_features_deprecated.cpp" )
create_single_source_cgal_program( "test_meshing_with_one_step.cpp" )
create_single_source_cgal_program( "test_meshing_with_one_step_with_features.cpp" )
create_single_source_cgal_program( "test_mesh_cell_base_3.cpp")
create_single_source_cgal_program( "test_min_edge_length.cpp")

foreach(target
    test_boost_has_xxx
    test_mesh_capsule_var_distance_bound
    test_implicit_multi_domain_to_labeling_function_wrapper
    test_criteria
    test_domain_with_polyline_features
    test_labeled_mesh_domain_3
    test_mesh_criteria_creation
    test_without_detect_features
    test_meshing_3D_image
    test_meshing_3D_image_with_features
    test_meshing_3D_gray_image
    test_meshing_implicit_function
    test_meshing_polyhedral_complex
    test_meshing_polyhedron
    test_meshing_polylines_only
    test_meshing_polyhedron_with_features
    test_meshing_verbose
    test_meshing_unit_tetrahedron
    test_meshing_with_default_edge_size
    test_meshing_determinism
    test_meshing_without_features_determinism
    test_mesh_3_issue_1554
    test_mesh_polyhedral_domain_with_features_deprecated
    test_mesh_cell_base_3
    test_meshing_with_one_step
    test_meshing_with_one_step_with_features
    test_min_edge_length
    test_min_size_criteria
    test_meshing_polyhedral_complex_with_manifold_and_min_size
    )
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
      test_meshing_determinism
      test_meshing_without_features_determinism
      test_mesh_3_issue_1554
      test_mesh_polyhedral_domain_with_features_deprecated
      test_mesh_cell_base_3
      test_min_edge_length
      test_min_size_criteria
      )
    if(TARGET ${target})
      target_link_libraries(${target} PUBLIC CGAL::TBB_support)
    endif()
  endforeach()

  if(CGAL_ENABLE_TESTING)
    set_property(TEST
      "execution   of  test_meshing_verbose"
      "execution   of  test_meshing_polyhedron_with_features"
      "execution   of  test_meshing_implicit_function"
      "execution   of  test_meshing_unit_tetrahedron"
      "execution   of  test_meshing_polyhedron"
      "execution   of  test_meshing_polyhedral_complex"
      "execution   of  test_mesh_capsule_var_distance_bound"
      "execution   of  test_mesh_3_issue_1554"
      "execution   of  test_mesh_polyhedral_domain_with_features_deprecated"
      "execution   of  test_mesh_cell_base_3"
      PROPERTY RUN_SERIAL 1)
    if(TARGET test_meshing_3D_image)
      set_property(TEST
        "execution   of  test_meshing_3D_image"
        "execution   of  test_meshing_3D_gray_image"
        PROPERTY RUN_SERIAL 1)
    endif()
  endif()
endif()

if(TARGET ITT::ITT)
  target_link_libraries(test_meshing_polyhedron_with_features PRIVATE ITT::ITT)
  target_compile_definitions(test_meshing_polyhedron_with_features PRIVATE CGAL_MESH_3_USE_INTEL_ITT)
  target_link_libraries(test_meshing_verbose PRIVATE ITT::ITT)
  target_compile_definitions(test_meshing_verbose PRIVATE CGAL_MESH_3_USE_INTEL_ITT)
endif()

if(CGAL_ENABLE_TESTING)
  set_tests_properties(
    "execution   of  test_meshing_polyhedron_with_features"
    "execution   of  test_meshing_verbose"
    PROPERTIES RESOURCE_LOCK Mesh_3_Tests_IO)
endif()
