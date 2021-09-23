cmake_minimum_required(VERSION 3.1...3.20)
project(Surface_mesh_topology_Tests)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)
if(CGAL_Qt5_FOUND)
  add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
endif()

# add_definitions(-DCGAL_TRACE_PATH_TESTS)
# add_definitions(-DCGAL_TRACE_CMAP_TOOLS)
# add_definitions("-Wall -Wextra")
# add_definitions("-D_GLIBCXX_DEBUG")

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

if(CGAL_Qt5_FOUND)
  target_link_libraries(fundamental_group_of_the_circle PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(fundamental_group_of_the_torus PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(path_with_rle_deformation_tests PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(test_homotopy PUBLIC CGAL::CGAL_Qt5)
endif()
