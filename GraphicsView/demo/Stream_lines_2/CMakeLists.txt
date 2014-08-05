# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Stream_lines_2_demo)

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )

find_package(CGAL COMPONENTS Qt4or5 REQUIRED)

include(${CGAL_USE_FILE})
	
find_package(Qt4or5)

set(CMAKE_AUTOMOC OFF)


include_directories (BEFORE ../../include)


if ( CGAL_FOUND AND QT${QT_VERSION_USED}_FOUND )

#--------------------------------
# UI files (Qt Designer files)
qt_wrap_ui( DT_UI_FILES Stream_lines_2.ui )

# qrc files (resources files, that contain icons, at least)
qt_add_resources ( DT_RESOURCE_FILES ./Stream_lines_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt_generate_moc( Stream_lines_2.cpp Stream_lines_2.moc )

# The executable itself.
add_executable  ( Stream_lines_2 Stream_lines_2.cpp Stream_lines_2.moc ${DT_UI_FILES} ${DT_RESOURCE_FILES} )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Stream_lines_2 )

# Link with Qt libraries
target_link_libraries( Stream_lines_2 ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries( Stream_lines_2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

else()

  message(STATUS "NOTICE: This demo requires CGAL and ${QT_VERSION}, and will not be compiled.")

endif()
