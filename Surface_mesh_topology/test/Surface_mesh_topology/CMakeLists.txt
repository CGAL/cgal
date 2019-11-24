project(Surface_mesh_topology_Tests)

cmake_minimum_required(VERSION 3.1)

find_package(CGAL COMPONENTS Qt5)
if(CGAL_Qt5_FOUND)
    add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
endif()

if ( NOT CGAL_FOUND )

  message(FATAL_ERROR "This program requires the CGAL library, and will not be compiled.")

endif()

set(HEADER_FILES
  ../../include/CGAL/Curves_on_surface_topology.h
  ../../include/CGAL/draw_face_graph_with_paths.h
  ../../include/CGAL/Path_on_surface.h
  ../../include/CGAL/Polygonal_schema.h
  ../../include/CGAL/Polygonal_schema_min_items.h
  ../../include/CGAL/Surface_mesh_topology/internal/Minimal_quadrangulation.h
  ../../include/CGAL/Surface_mesh_topology/internal/Path_generators.h
  ../../include/CGAL/Surface_mesh_topology/internal/Path_on_surface_with_rle.h
  )

set(SOURCE_FILES
  fundamental_group_of_the_circle.cpp
  fundamental_group_of_the_torus.cpp
  homotopy_big_cylinder.cpp
  homotopy_double_torus_with_holes.cpp
  homotopy_rond_point_saucisse.cpp
  homotopy_small_cylinder.cpp
  path_tests.cpp
  path_with_rle_deformation_tests.cpp
  test_homotopy.cpp
  test_homotopy_with_polygonal_schema.cpp
  tests_path_on_surface.cpp
  )

foreach(cppfile ${SOURCE_FILES})
  create_single_source_cgal_program("${cppfile}" ${HEADER_FILES})
endforeach()

if(CGAL_Qt5_FOUND)
  target_link_libraries(fundamental_group_of_the_circle PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(fundamental_group_of_the_torus PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(path_with_rle_deformation_tests PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(test_homotopy PUBLIC CGAL::CGAL_Qt5)
endif()
