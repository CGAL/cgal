# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (GraphicsView)

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL REQUIRED Qt4)

include(${CGAL_USE_FILE})

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )


find_package(Qt4       REQUIRED)
include(${QT_USE_FILE})
include_directories (BEFORE ../../include)
include_directories (BEFORE ./include)
link_directories(../../src/CGALQt4)

#--------------------------------
# The "constrained Delaunay" demo: Constrained_Delaunay_triangulation_2
#--------------------------------
# UI files (Qt Designer files)
qt4_wrap_ui( CDT_UI_FILES Constrained_Delaunay_triangulation_2.ui )

# qrc files (resources files, that contain icons, at least)
qt4_add_resources ( CDT_RESOURCE_FILES ./Constrained_Delaunay_triangulation_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt4_generate_moc( Constrained_Delaunay_triangulation_2.cpp Constrained_Delaunay_triangulation_2.moc )

# The executable itself.
add_executable  ( Constrained_Delaunay_triangulation_2 Constrained_Delaunay_triangulation_2.cpp Constrained_Delaunay_triangulation_2.moc ${CDT_UI_FILES} ${CDT_RESOURCE_FILES} )

# Link with Qt libraries
target_link_libraries( Constrained_Delaunay_triangulation_2 ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries( Constrained_Delaunay_triangulation_2 ${CGAL_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES})

#--------------------------------
# The "Delaunay" demo: Delaunay_triangulation_2
#--------------------------------
# UI files (Qt Designer files)
qt4_wrap_ui( DT_UI_FILES Delaunay_triangulation_2.ui )

# qrc files (resources files, that contain icons, at least)
qt4_add_resources ( DT_RESOURCE_FILES ./Delaunay_triangulation_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt4_generate_moc( Delaunay_triangulation_2.cpp Delaunay_triangulation_2.moc )

# The executable itself.
add_executable  ( Delaunay_triangulation_2 Delaunay_triangulation_2.cpp Delaunay_triangulation_2.moc ${DT_UI_FILES} ${DT_RESOURCE_FILES} )

# Link with Qt libraries
target_link_libraries( Delaunay_triangulation_2 ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries( Delaunay_triangulation_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
