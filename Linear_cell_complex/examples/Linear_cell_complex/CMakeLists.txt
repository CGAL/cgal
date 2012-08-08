# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Linear_cell_complex_examples )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL QUIET COMPONENTS Core Qt4)

# For Gprof.
# ADD_DEFINITIONS("-pg")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

# To use valgrind, we must disable rounding math ckeck.
# add_definition(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )
  include_directories(BEFORE ../../include)

  create_single_source_cgal_program( "linear_cell_complex_3.cpp" )
  create_single_source_cgal_program( "linear_cell_complex_4.cpp" )
  create_single_source_cgal_program(
                  "linear_cell_complex_3_with_colored_vertices.cpp" )
  create_single_source_cgal_program(
                  "linear_cell_complex_3_with_mypoint.cpp" )

  create_single_source_cgal_program("plane_graph_to_lcc_2.cpp")

# If you want to visualize a linear cell complex, there are 2 viewers 
# based on qt and vtk. Just uncomment the corresponding line
  #  include("CMakeLCCViewerQt.inc")
  #  include("CMakeLCCViewerVtk.inc")

  add_executable(linear_cell_complex_3_triangulation
                          linear_cell_complex_3_triangulation.cpp)
  target_link_libraries(linear_cell_complex_3_triangulation ${CGAL_LIBRARIES}
                                   ${CGAL_3RD_PARTY_LIBRARIES}
                                   ${MAP_VIEWER_LIBRARIES})

  add_executable(voronoi_2 voronoi_2.cpp)
  target_link_libraries(voronoi_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}
                                  ${MAP_VIEWER_LIBRARIES})

  add_executable(voronoi_3 voronoi_3.cpp)
  target_link_libraries(voronoi_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}
                                  ${MAP_VIEWER_LIBRARIES})

else()
  
  message(STATUS "This program requires the CGAL library, "
                  "and will not be compiled.")
  
endif()

