# This is the CMake script for compiling a CGAL application.

project( Apollonius_graph_2 ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL REQUIRED Qt3 )
include( ${CGAL_USE_FILE} )

find_package(Qt3 REQUIRED)

include_directories (BEFORE ./include)

# use the Qt MOC preprocessor on classes that derives from QObject
include( Qt3Macros )
qt3_generate_moc( apollonius_graph_2.cpp apollonius_graph_2.moc )
qt3_generate_moc( qt_file_toolbar.h qt_file_toolbar.moc )
qt3_generate_moc( qt_layers_toolbar.h qt_layers_toolbar.moc )
qt3_generate_moc( edit_vertex_layer.h edit_vertex_layer.moc )


# The executable itself.
add_executable  ( apollonius_graph_2
		  apollonius_graph_2.cpp apollonius_graph_2.moc 
		  qt_file_toolbar.moc qt_layers_toolbar.moc edit_vertex_layer.moc )


target_link_libraries( apollonius_graph_2  ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
