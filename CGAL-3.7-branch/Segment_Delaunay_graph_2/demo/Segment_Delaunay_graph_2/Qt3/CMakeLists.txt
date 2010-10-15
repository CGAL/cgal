# This is the CMake script for compiling a CGAL application.

project( Segment_Delaunay_graph_2 ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL COMPONENTS Qt3 Core )
include( ${CGAL_USE_FILE} )

find_package(Qt3-patched)
# FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
# that it can be used together with FindQt4: all its variables are prefixed
# by "QT3_" instead of "QT_".

if ( CGAL_FOUND AND CGAL_Qt3_FOUND AND QT3_FOUND )

  include_directories (BEFORE ./include)

  # use the Qt MOC preprocessor on classes that derives from QObject
  include( Qt3Macros-patched )

  qt3_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/qt_file_toolbar.h" qt_file_toolbar.moc )
  qt3_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/qt_layers_toolbar.h" qt_layers_toolbar.moc )
  qt3_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/my_window.h" my_window.moc )


  # The executable itself.
  add_executable  ( polygon_delaunay_graph_2
                    polygon_delaunay_graph_2.cpp  
                    qt_file_toolbar.moc  qt_layers_toolbar.moc  my_window.moc )

  add_executable  ( segment_delaunay_graph_2
                    segment_delaunay_graph_2.cpp
                    qt_file_toolbar.moc  qt_layers_toolbar.moc  my_window.moc )

  target_link_libraries( polygon_delaunay_graph_2  ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  target_link_libraries( segment_delaunay_graph_2  ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS polygon_delaunay_graph_2 )
  add_to_cached_list( CGAL_EXECUTABLE_TARGETS segment_delaunay_graph_2 )
else()

  message(STATUS "NOTICE: This demo requires CGAL,, Core, and Qt3, and will not be compiled.")

endif()
