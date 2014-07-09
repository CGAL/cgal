cmake_minimum_required(VERSION 2.4.5)

project (Periodic_2_Triangulation_2_demo)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
endif()


set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )


#New for Qt5 version !
option (USE_QT5 "Use Qt5 version instead of Qt4" ON)

if(USE_QT5)

	find_package(CGAL COMPONENTS Qt5 REQUIRED)

	include(${CGAL_USE_FILE})
	
	find_package(Qt5)

        set(CMAKE_AUTOMOC OFF)

else (USE_QT5)

	find_package(CGAL COMPONENTS Qt4 REQUIRED)

	include(${CGAL_USE_FILE})

	find_package(Qt4)

	if(CGAL_Qt4_FOUND AND QT4_FOUND)
		include(${QT_USE_FILE})
		message("Qt4 found")
		set(QT4 TRUE)
	endif()

endif (USE_QT5)

include(QtChoice)


include_directories (BEFORE ../../../include)
include_directories (BEFORE ./include)

if ( QT4 OR QT5 )
#--------------------------------
# The "2D Periodic triangulation" demo: Periodic_2_triangulation_2
#--------------------------------
# UI files (Qt Designer files)
qt_wrap_ui( DT_UI_FILES Periodic_2_triangulation_2.ui )

# qrc files (resources files, that contain icons, at least)
qt_add_resources ( DT_RESOURCE_FILES ./Periodic_2_triangulation_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt_generate_moc( Periodic_2_Delaunay_triangulation_2.cpp Periodic_2_triangulation_2.moc )

# find header files for projects that can show them
file(GLOB headers "*.h")
file(GLOB QT_headers "include/CGAL/Qt/*.h")
file(GLOB P2T2_headers "../../../include/CGAL/*.h")
SOURCE_GROUP("P2T2" FILES ${P2T2_headers})
SOURCE_GROUP("QT" FILES ${QT_headers})

# The executable itself.
add_executable  ( Periodic_2_Delaunay_triangulation_2 
  Periodic_2_Delaunay_triangulation_2.cpp
  Periodic_2_triangulation_2.moc 
  ${DT_UI_FILES} ${DT_RESOURCE_FILES} 
  ${headers} ${QT_headers} ${P2T2_headers})

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Periodic_2_Delaunay_triangulation_2 )

# Link with Qt libraries
target_link_libraries( Periodic_2_Delaunay_triangulation_2 ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries( Periodic_2_Delaunay_triangulation_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt4 or Qt5, and will not be compiled.")

endif()
