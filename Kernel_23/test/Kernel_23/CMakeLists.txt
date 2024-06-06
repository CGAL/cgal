cmake_minimum_required(VERSION 3.12...3.29)
project(Kernel_23_Tests)

find_package(CGAL REQUIRED COMPONENTS Core)

include_directories(BEFORE "include")

create_single_source_cgal_program("Cartesian.cpp")
create_single_source_cgal_program("determinant_77.cpp")
create_single_source_cgal_program("Dimension.cpp")
create_single_source_cgal_program("Exact_predicates_exact_constructions_kernel.cpp")
create_single_source_cgal_program("Filtered_cartesian.cpp")
create_single_source_cgal_program("Filtered_homogeneous.cpp")
create_single_source_cgal_program("Homogeneous.cpp")
create_single_source_cgal_program("issue_129.cpp")
create_single_source_cgal_program("issue_3301.cpp")
create_single_source_cgal_program("issue_8140.cpp")
create_single_source_cgal_program("Kernel_checker.cpp")
create_single_source_cgal_program("Lazy_kernel.cpp")
create_single_source_cgal_program("origin_3.cpp")
create_single_source_cgal_program("overload_bug.cpp")
create_single_source_cgal_program("rank.cpp")
create_single_source_cgal_program("Simple_cartesian.cpp")
create_single_source_cgal_program("Simple_homogeneous.cpp")
create_single_source_cgal_program("test_all_linear_intersections.cpp")
create_single_source_cgal_program("test_approximate_dihedral_angle_3.cpp")
create_single_source_cgal_program("test_bbox.cpp")
create_single_source_cgal_program("test_converter.cpp")
create_single_source_cgal_program("test_Has_conversion.cpp")
create_single_source_cgal_program("test_hash_functions.cpp")
create_single_source_cgal_program("test_kernel__.cpp")
create_single_source_cgal_program("test_projection_traits.cpp")
create_single_source_cgal_program("test_Projection_traits_xy_3_Intersect_2.cpp")

if(CGAL_KERNEL_23_TEST_RT_FT_PREDICATE_FLAGS)
  # Templated operators:
  # - create a lot of possible combinations, which is expensive to test
  # - create issues because some combinations might be RT-sufficient whereas others will require FT
  #
  # add_definitions(-DCGAL_KERNEL_23_TEST_RT_FT_PREDICATES_TEST_PREDICATES_WITH_TEMPLATED_OPERATORS)

  create_single_source_cgal_program("atomic_compilation_test.cpp")
  target_precompile_headers(atomic_compilation_test PUBLIC [["atomic_RT_FT_predicate_headers.h"]])

  create_single_source_cgal_program("test_RT_or_FT_predicates.cpp")
  target_compile_definitions(test_RT_or_FT_predicates PRIVATE
    "CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
    "CMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}")
endif()
