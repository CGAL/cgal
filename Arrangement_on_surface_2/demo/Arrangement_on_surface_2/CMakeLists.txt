# This is the CMake script for compiling a CGAL application.

project( Arrangement_2 ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL REQUIRED Qt3 )
include( ${CGAL_USE_FILE} )

find_package(Qt3 REQUIRED)

# use the Qt MOC preprocessor on classes that derives from QObject
include( Qt3Macros )
qt3_generate_moc( arrangement_2.h arrangement_2.moc )
qt3_generate_moc( forms.h forms.moc )


# The executable itself.
add_executable  ( arrangement_2 
		  arrangement_2.cpp arrangement_2.moc
                  forms.cpp forms.moc
		  demo_tab.cpp
		  MyWindow_files.cpp
		  MyWindow_operations.cpp
		  at_layer.cpp
		  )


