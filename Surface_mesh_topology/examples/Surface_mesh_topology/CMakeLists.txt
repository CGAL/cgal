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
  edgewidth_lcc.cpp
  edgewidth_surface_mesh.cpp
  facewidth.cpp
  map_2_constructor.cpp
  open_path_homotopy.cpp
  path_homotopy.cpp
  path_homotopy_double_torus.cpp
  path_homotopy_torus.cpp
  path_homotopy_with_sm_and_polyhedron.cpp
  path_homotopy_with_symbols_2.cpp
  path_homotopy_with_symbols.cpp
  shortest_noncontractible_cycle_2.cpp
  shortest_noncontractible_cycle.cpp
  unsew_edgewidth_repeatedly.cpp
  )

set(HEADER_FILES
  ../../include/CGAL/Curves_on_surface_topology.h
  ../../include/CGAL/Polygonal_schema_fwd.h
  ../../include/CGAL/draw_face_graph_with_paths.h
  ../../include/CGAL/Polygonal_schema_min_items.h
  ../../include/CGAL/Polygonal_schema.h
  ../../include/CGAL/Path_on_surface.h
  ../../include/CGAL/Face_graph_wrapper.h
  ../../include/CGAL/Surface_mesh_topology/internal/Minimal_quadrangulation.h
  ../../include/CGAL/Surface_mesh_topology/internal/Generic_map_selector.h
  ../../include/CGAL/Surface_mesh_topology/internal/Path_on_surface_with_rle.h
  ../../include/CGAL/Surface_mesh_topology/internal/Facewidth.h
  ../../include/CGAL/Surface_mesh_topology/internal/Shortest_noncontractible_cycle.h
  ../../include/CGAL/Surface_mesh_topology/internal/Path_generators.h
  ../../include/CGAL/Surface_mesh_topology/internal/Functors_for_face_graph_wrapper.h
  ../../include/CGAL/Surface_mesh_topology/internal/Iterators_for_face_graph_wrapper.h
  )

foreach(cppfile ${SOURCE_FILES})
  create_single_source_cgal_program("${cppfile}" "${HEADER_FILES}")
endforeach()

if(CGAL_Qt5_FOUND)
  target_link_libraries(edgewidth_lcc PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(edgewidth_surface_mesh PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(facewidth PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(open_path_homotopy PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(path_homotopy PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(path_homotopy_double_torus PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(path_homotopy_torus PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(path_homotopy_with_sm_and_polyhedron PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(shortest_noncontractible_cycle_2 PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(shortest_noncontractible_cycle PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(unsew_edgewidth_repeatedly PUBLIC CGAL::CGAL_Qt5)
endif()
