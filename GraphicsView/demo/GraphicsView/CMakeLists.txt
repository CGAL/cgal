# This is the CMake script for compiling a CGAL application.

project (GraphicsView)

cmake_minimum_required(VERSION 2.4)
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL REQUIRED)

include_directories (BEFORE ../../include)
include_directories (BEFORE ./include)
link_directories(../../src/CGALQt4)

include(${CGAL_USE_FILE})

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )


find_package(Qt4       REQUIRED)
include(${QT_USE_FILE})
include_directories (BEFORE ../../include)
link_directories(../../src/CGALQt4)

#----------------------------------------
# minimal graphics view executable: min 
#----------------------------------------
add_executable  ( min min.cpp )
if ( NOT AUTO_LINK_ENABLED )    
  # If auto-link is not enabled (on Linux for example), link with libCGAL explicitly.
  target_link_libraries( min ${QT_LIBRARIES} CGAL CGAL-Qt4 )
endif()

