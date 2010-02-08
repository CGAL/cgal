# This is the CMake script for compiling a CGAL application.

project( Apollonius_graph_2 ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL COMPONENTS Qt3 )
include( ${CGAL_USE_FILE} )

find_package(Qt3-patched)
# FindQt3-patched.cmake is FindQt3.cmake patched by CGAL developers, so
# that it can be used together with FindQt4: all its variables are prefixed
# by "QT3_" instead of "QT_".

if ( CGAL_FOUND AND QT3_FOUND AND CGAL_Qt3_FOUND )

  include_directories (BEFORE ./include)
  
  # use the Qt MOC preprocessor on classes that derives from QObject
  include( Qt3Macros-patched )
  qt3_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/apollonius_graph_2.cpp" apollonius_graph_2.moc )
  qt3_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/qt_file_toolbar.h" qt_file_toolbar.moc )
  qt3_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/qt_layers_toolbar.h" qt_layers_toolbar.moc )
  qt3_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/edit_vertex_layer.h" edit_vertex_layer.moc )
  
  
  # The executable itself.
  add_executable  ( apollonius_graph_2
  		  apollonius_graph_2.cpp apollonius_graph_2.moc 
  		  qt_file_toolbar.moc qt_layers_toolbar.moc edit_vertex_layer.moc )
  
  add_to_cached_list( CGAL_EXECUTABLE_TARGETS apollonius_graph_2 )
  
  target_link_libraries( apollonius_graph_2  ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt3, and will not be compiled.")

endif()
