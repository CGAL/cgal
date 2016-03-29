# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Linear_cell_complex_examples )

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

# If you want to visualize a linear cell complex, you can use the following viewer
# based on qt. Just uncomment the following two lines, plus the lines qt5_use_modules below

# find_package(CGAL COMPONENTS Qt5)
# include("CMakeLCCViewerQt.inc")

# If you don't want to visualize, use the following line (otherwise comment it)
find_package(CGAL QUIET)

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
  create_single_source_cgal_program("linear_cell_complex_3_attributes_management.cpp")

  add_executable(voronoi_2 voronoi_2.cpp)
  target_link_libraries(voronoi_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}
                                  ${MAP_VIEWER_LIBRARIES})
  # qt5_use_modules(voronoi_2 ${MAP_VIEWER_MODULES})

  add_executable(voronoi_3 voronoi_3.cpp)
  target_link_libraries(voronoi_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}
                                  ${MAP_VIEWER_LIBRARIES})
  # qt5_use_modules(voronoi_3 ${MAP_VIEWER_MODULES})

else()

  message(STATUS "This program requires the CGAL library, "
                  "and will not be compiled.")

endif()

