# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Linear_cell_complex_Examples )

cmake_minimum_required(VERSION 3.1)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

find_package(CGAL COMPONENTS Qt5)

# For Gprof.
# ADD_DEFINITIONS("-pg")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

# To use valgrind, we must disable rounding math ckeck.
# add_definition(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)

if ( CGAL_FOUND )
  include(${CGAL_USE_FILE})

  include( CGAL_CreateSingleSourceCGALProgram )
  include_directories(BEFORE ../../include)

  create_single_source_cgal_program("linear_cell_complex_3.cpp")
  create_single_source_cgal_program("linear_cell_complex_4.cpp")
  create_single_source_cgal_program("linear_cell_complex_3_with_colored_vertices.cpp")
  create_single_source_cgal_program("linear_cell_complex_3_with_mypoint.cpp")
  create_single_source_cgal_program("plane_graph_to_lcc_2.cpp")
  create_single_source_cgal_program("linear_cell_complex_3_attributes_management.cpp")
  create_single_source_cgal_program("voronoi_2.cpp")
  create_single_source_cgal_program("voronoi_3.cpp")
  create_single_source_cgal_program("gmap_linear_cell_complex_3.cpp")
  create_single_source_cgal_program("linear_cell_complex_3_operations.cpp")  

  create_single_source_cgal_program("draw_linear_cell_complex.cpp")                           
  if(TARGET CGAL::CGAL_Qt5)
    find_package(Qt5        REQUIRED COMPONENTS OpenGL Xml)
    find_package(QGLViewer  REQUIRED)
    find_package(OpenGL     REQUIRED)
    add_definitions(${QT_DEFINITIONS})
    add_definitions(-DQT_NO_KEYWORDS)
    include_directories( ${QGLVIEWER_INCLUDE_DIR} )
    add_definitions(${QGLVIEWER_DEFINITIONS})

    target_link_libraries(draw_linear_cell_complex CGAL::CGAL_Qt5 Qt5::Xml
                                                      ${QT_LIBRARIES}
                                                      ${QGLVIEWER_LIBRARIES}
                                                      ${OPENGL_gl_LIBRARY})
    target_compile_definitions(draw_linear_cell_complex PRIVATE CGAL_USE_BASIC_VIEWER)
  endif()

else()

  message(STATUS "This program requires the CGAL library, "
                  "and will not be compiled.")

endif()

