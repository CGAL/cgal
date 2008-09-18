# This is the CMake script for compiling a CGAL application.

project( Boolean_operations_2 ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL REQUIRED Qt3 )

include( ${CGAL_USE_FILE} )


# use the Qt MOC preprocessor on classes that derives from QObject
include( Qt3Macros )
qt3_generate_moc( boolean_operations_2_toolbar.h boolean_operations_2_toolbar.moc )
qt3_generate_moc( boolean_operations_2.h boolean_operations_2.moc )


# The executable itself.
add_executable  ( boolean_operations_2 
		  boolean_operations_2.cpp boolean_operations_2.moc
                  boolean_operations_2_toolbar.cpp boolean_operations_2_toolbar.moc
                  Qt_widget_locate_layer.cpp
		  )


