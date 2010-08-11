# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Bounding_volumes)

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL COMPONENTS Qt4)

include(${CGAL_USE_FILE})

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )


find_package(Qt4)
include_directories (BEFORE ../../include)
include_directories (BEFORE ./include)
include_directories (BEFORE ../../../Matrix_search/include)


if ( CGAL_FOUND AND CGAL_Qt4_FOUND AND QT4_FOUND )

include(${QT_USE_FILE})

#----------------------------------------------
# The "Bounding volumes" demo: Bounding_volumes
#----------------------------------------------
# UI files (Qt Designer files)
qt4_wrap_ui( DT_UI_FILES Bounding_volumes.ui )

# qrc files (resources files, that contain icons, at least)
qt4_add_resources ( DT_RESOURCE_FILES ./Bounding_volumes.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt4_generate_moc( "Bounding_volumes.cpp" "${CMAKE_CURRENT_BINARY_DIR}/Bounding_volumes.moc" )

# The executable itself.
add_executable  ( Bounding_volumes Bounding_volumes.cpp Bounding_volumes.moc ${DT_UI_FILES} ${DT_RESOURCE_FILES} )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS Bounding_volumes )

# Link with Qt libraries
target_link_libraries( Bounding_volumes ${QT_LIBRARIES} )
# Link with CGAL
target_link_libraries( Bounding_volumes ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt4, and will not be compiled.")

endif()
