# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Segment_voronoi_2)

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL COMPONENTS Qt4 Core)

include(${CGAL_USE_FILE})

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )


find_package(Qt4)
include_directories (BEFORE ../../include)
include_directories (BEFORE ./include)


if ( CGAL_FOUND AND CGAL_Qt4_FOUND AND QT4_FOUND )

include(${QT_USE_FILE})

#--------------------------------
# The "Segment Voronoi" demo: Segment_voronoi_2
#--------------------------------
# UI files (Qt Designer files)
qt4_wrap_ui( CDT_UI_FILES Segment_voronoi_2.ui )

# qrc files (resources files, that contain icons, at least)
qt4_add_resources ( CDT_RESOURCE_FILES ./Segment_voronoi_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt4_generate_moc( Segment_voronoi_2.cpp Segment_voronoi_2.moc )

# The executable itself.
add_executable  ( Segment_voronoi_2 Segment_voronoi_2.cpp Segment_voronoi_2.moc ${CDT_UI_FILES} ${CDT_RESOURCE_FILES} )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Segment_voronoi_2 )

# Link with Qt libraries
target_link_libraries( Segment_voronoi_2 ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries( Segment_voronoi_2 ${CGAL_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES})


else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt4, and will not be compiled.")

endif()
