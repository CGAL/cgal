# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

# cmake ../ -DCMAKE_BUILD_TYPE=Debug 
# ou
# cmake ../ -DCMAKE_BUILD_TYPE=Release

project( Map_3_examples )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall")
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
endif()

# ADD_DEFINITIONS("-pg")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

# Pour le problème de valgrind avec CGAL
# add_definition(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )
  include_directories(BEFORE ../../include)
  include_directories(BEFORE ../../../Combinatorial_map/include)

  create_single_source_cgal_program( "linear_cell_complex_3.cpp" )
  create_single_source_cgal_program( "linear_cell_complex_4.cpp" )
  create_single_source_cgal_program( "linear_cell_complex_3_with_colored_vertices.cpp" )
  create_single_source_cgal_program( "map_3_subdivision.cpp" )
  create_single_source_cgal_program( "plane_graph_to_map_2.cpp" )
  create_single_source_cgal_program( "map_3_iterators.cpp" )
  create_single_source_cgal_program( "exemple_incremental_builder.cpp" )
  create_single_source_cgal_program( "polyhedron_clear.cpp" )

# If you want to visualize a map, there are 2 viewers based on qt and vtk
  #include_directories("../../include/CGAL/Combinatorial_map_viewers/" )
  #include("../../include/CGAL/Combinatorial_map_viewers/CMakeMapViewerQt.inc")
  #include("../../include/CGAL/Combinatorial_map_viewers/CMakeMapViewerVtk.inc")

  add_executable(voronoi_3 voronoi_3.cpp)
  target_link_libraries(voronoi_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
  # target_link_libraries(voronoi_3 ${MAP_VIEWER_LIBRARIES_QT})
  # target_link_libraries(voronoi_3 ${MAP_VIEWER_LIBRARIES_VTK}) 
else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

