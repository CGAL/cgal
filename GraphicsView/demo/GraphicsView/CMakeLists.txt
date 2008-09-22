# This is the CMake script for compiling a CGAL application.

project (GraphicsView)

cmake_minimum_required(VERSION 2.4)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL REQUIRED Qt4)

include(${CGAL_USE_FILE})

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )

find_package(Qt4       REQUIRED)
include(${QT_USE_FILE})


add_executable  ( min min.cpp )


# Link with Qt libraries
target_link_libraries( min ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries( min ${CGAL_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES})