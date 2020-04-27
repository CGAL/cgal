cmake_minimum_required(VERSION 3.1...3.15)
project(Surface_mesh_topology_Examples)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

find_package(CGAL QUIET COMPONENTS Qt5)
if(CGAL_Qt5_FOUND)
    add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
endif()

if ( NOT CGAL_FOUND )
  message(FATAL_ERROR "This program requires the CGAL library, and will not be compiled.")
endif()

# add_definitions(-DCGAL_TRACE_PATH_TESTS)
# add_definitions(-DCGAL_TRACE_CMAP_TOOLS)
# add_definitions("-Wall -Wextra")
# add_definitions("-D_GLIBCXX_DEBUG")

set(SOURCE_FILES
  basic_example_surface_mesh_topology.cpp
  basic_example_torus.cpp
  map_2_constructor.cpp
  open_path_homotopy.cpp
  path_homotopy.cpp
  path_homotopy_with_symbols.cpp
  path_homotopy_with_symbols_2.cpp
  surface_mesh_topology_with_sm_and_polyhedron.cpp
  shortest_noncontractible_cycle_through_a_vertex.cpp
  shortest_noncontractible_cycle_using_BFS.cpp
  unsew_edgewidth_repeatedly.cpp
  edgewidth_surface_mesh.cpp
  facewidth_on_unweighted_map.cpp
  )

foreach(cppfile ${SOURCE_FILES})
  create_single_source_cgal_program("${cppfile}")
endforeach()

if(CGAL_Qt5_FOUND)
  target_link_libraries(basic_example_surface_mesh_topology PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(basic_example_torus PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(open_path_homotopy PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(path_homotopy PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(surface_mesh_topology_with_sm_and_polyhedron PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(unsew_edgewidth_repeatedly PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(shortest_noncontractible_cycle_through_a_vertex PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(shortest_noncontractible_cycle_using_BFS PUBLIC CGAL::CGAL_Qt5)
endif()
