# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Linear_cell_complex_Examples )

cmake_minimum_required(VERSION 3.1)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

# If you want to visualize a linear cell complex, you can use the following viewer
# based on qt. Just uncomment the following two lines, plus the lines qt5_use_modules below

# find_package(CGAL COMPONENTS Qt5)
# include("CMakeBasicViewerQt.inc")

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
                                  ${BASIC_VIEWER_LIBRARIES})
  if (USE_BASIC_VIEWER)
    qt5_use_modules(voronoi_2 ${BASIC_VIEWER_MODULES})
  endif(USE_BASIC_VIEWER)

  add_executable(voronoi_3 voronoi_3.cpp)
  target_link_libraries(voronoi_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}
                                  ${BASIC_VIEWER_LIBRARIES})

  if (USE_BASIC_VIEWER)
    qt5_use_modules(voronoi_3 ${BASIC_VIEWER_MODULES})
  endif(USE_BASIC_VIEWER)
    
  create_single_source_cgal_program( "gmap_linear_cell_complex_3.cpp" )

  add_executable(linear_cell_complex_3_operations linear_cell_complex_3_operations.cpp)
  target_link_libraries(linear_cell_complex_3_operations ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES}
                                  ${BASIC_VIEWER_LIBRARIES})
  if (USE_BASIC_VIEWER)
    qt5_use_modules(linear_cell_complex_3_operations ${BASIC_VIEWER_MODULES})
  endif(USE_BASIC_VIEWER)
  
else()

  message(STATUS "This program requires the CGAL library, "
                  "and will not be compiled.")

endif()

