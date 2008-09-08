# This is the CMake script for compiling a CGAL application.

project( Alpha_shapes_2 ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL REQUIRED Qt3 )
include( ${CGAL_USE_FILE} )

find_package(Qt3 REQUIRED)

# use the Qt MOC preprocessor on classes that derives from QObject
include( Qt3Macros )
qt3_generate_moc( Qt_widget_move_point.h Qt_widget_move_point.moc )
qt3_generate_moc( Qt_widget_toolbar.h Qt_widget_toolbar.moc )
qt3_generate_moc( Qt_widget_toolbar_layers.h Qt_widget_toolbar_layers.moc )
qt3_generate_moc( alpha_shapes_2.cpp alpha_shapes_2.moc )


# The executable itself.
add_executable  ( alpha_shapes_2 
		  alpha_shapes_2.cpp alpha_shapes_2.moc
		  Qt_widget_move_point.cpp Qt_widget_move_point.moc
		  Qt_widget_toolbar.cpp Qt_widget_toolbar.moc 
		  Qt_widget_toolbar_layers.cpp Qt_widget_toolbar_layers.moc )


