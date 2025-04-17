cmake_minimum_required(VERSION 3.12...3.31)
project(Surface_mesh_topology_Tests)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

# add_compile_definitions(CGAL_TRACE_PATH_TESTS)
# add_compile_definitions(CGAL_TRACE_CMAP_TOOLS)
# add_definitions("-Wall -Wextra")
# add_compile_definitions(_GLIBCXX_DEBUG)

set(SOURCE_FILES
  fundamental_group_of_the_circle.cpp
  fundamental_group_of_the_torus.cpp
  homotopy_big_cylinder.cpp
  homotopy_double_torus_with_holes.cpp
  homotopy_rond_point_saucisse.cpp
  homotopy_small_cylinder.cpp
  path_tests.cpp
  path_with_rle_deformation_tests.cpp
  shortest_noncontractible_cycle_tests.cpp
  simplicity_cylinder.cpp
  simplicity_double_torus.cpp
  simplicity_double_torus_with_holes.cpp
  simplicity_fundamental_polygon.cpp
  simplicity_torus.cpp
  test_homotopy.cpp
  test_homotopy_with_polygonal_schema.cpp
  test_shortest_cycle_non_contractible.cpp
  tests_path_on_surface.cpp
  )

foreach(cppfile ${SOURCE_FILES})
  create_single_source_cgal_program("${cppfile}")
endforeach()

if(CGAL_Qt6_FOUND)
  target_link_libraries(fundamental_group_of_the_circle PRIVATE CGAL::CGAL_Basic_viewer_Qt)
  target_link_libraries(fundamental_group_of_the_torus PRIVATE CGAL::CGAL_Basic_viewer_Qt)
  target_link_libraries(path_with_rle_deformation_tests PRIVATE CGAL::CGAL_Basic_viewer_Qt)
  target_link_libraries(test_homotopy PRIVATE CGAL::CGAL_Basic_viewer_Qt)
endif()
