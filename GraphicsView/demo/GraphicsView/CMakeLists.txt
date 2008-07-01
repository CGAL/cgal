# This is the CMake script for compiling a CGAL application.

project (GraphicsView)

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

#----------------------------------------
# minimal graphics view executable: min 
#----------------------------------------
add_executable  ( min EXCLUDE_FROM_ALL min.cpp )
qt4_generate_moc( min.cpp min.moc )
target_link_libraries( min ${QT_LIBRARIES} )

