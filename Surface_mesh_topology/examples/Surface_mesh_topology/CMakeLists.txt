project(Surface_mesh_topology_Examples)

cmake_minimum_required(VERSION 3.1)

find_package(CGAL QUIET COMPONENTS Qt5)
if(CGAL_Qt5_FOUND)
    add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
endif()

# add_definitions(-DCGAL_TRACE_PATH_TESTS)
# add_definitions(-DCGAL_TRACE_CMAP_TOOLS)

if ( NOT CGAL_FOUND )

  message(FATAL_ERROR "This program requires the CGAL library, and will not be compiled.")

endif()

include_directories(../../include)

set(SOURCE_FILES
  # basic_example_surface_mesh_topology.cpp
  # basic_example_torus.cpp
  # map_2_constructor.cpp
  # open_path_homotopy.cpp
  path_homotopy.cpp
  # path_homotopy_with_symbols_2.cpp
  # path_homotopy_with_symbols.cpp
  # surface_mesh_topology_with_sm_and_polyhedron.cpp
   )

foreach(cppfile ${SOURCE_FILES})
  create_single_source_cgal_program("${cppfile}")
endforeach()

if(CGAL_Qt5_FOUND)
  # target_link_libraries(basic_example_surface_mesh_topology PUBLIC CGAL::CGAL_Qt5)
  # target_link_libraries(basic_example_torus PUBLIC CGAL::CGAL_Qt5)
  # target_link_libraries(open_path_homotopy PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(path_homotopy PUBLIC CGAL::CGAL_Qt5)
  # target_link_libraries(surface_mesh_topology_with_sm_and_polyhedron PUBLIC CGAL::CGAL_Qt5)
endif()
