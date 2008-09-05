# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (GraphicsView)

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL REQUIRED)

include(${CGAL_USE_FILE})

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )


find_package(Qt4       REQUIRED)
include(${QT_USE_FILE})
include_directories (BEFORE ../../include)
link_directories(../../src/CGALQt4)


#--------------------------------
# Demo: Polygon_2
#--------------------------------
# UI files (Qt Designer files)
qt4_wrap_ui( DT_UI_FILES Polygon_2.ui )

# qrc files (resources files, that contain icons, at least)
qt4_add_resources ( DT_RESOURCE_FILES ./Polygon_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt4_generate_moc( Polygon_2.cpp Polygon_2.moc )

# The executable itself.
add_executable  ( Polygon_2 Polygon_2.cpp ${DT_UI_FILES} ${DT_RESOURCE_FILES} )

# Link with Qt libraries
target_link_libraries( Polygon_2 ${QT_LIBRARIES} )
if ( NOT AUTO_LINK_ENABLED )    
  # If auto-link is not enabled (on Linux for example), link with libCGAL and libCGAL-Qt4 explicitly.
  target_link_libraries( Polygon_2 CGAL )
  target_link_libraries( Polygon_2 CGAL-Qt4 )
endif()
