# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

# ** DEVELOPMENT**  CMakeLists.txt.. 
# It requires the env var $CGAL_SVN with the location of a svn working copy
# And $CGAL_DIR with the location of a CGAL installation

project (Polyline_simplification_2)

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

if ( NOT CGAL_SVN )
  set( CGAL_SVN_env "$ENV{CGAL_SVN}" )
  if ( NOT "${CGAL_SVN_env}" STREQUAL "" )
    set( CGAL_SVN ${CGAL_SVN_env} CACHE FILEPATH "Location of the CGAL svn working copy" )
  else()  
    set( CGAL_SVN "NOTFOUND" CACHE FILEPATH "Location of the CGAL svn working copy" )
  endif()
endif()

include_directories (BEFORE ../../../../../trunk/GraphicsView/include)

# include_directories (BEFORE ${CGAL_SVN}/GraphicsView/include)
include_directories (BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories (BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/../../include )
include_directories (BEFORE ${CMAKE_CURRENT_BINARY_DIR})

if ( CGAL_FOUND AND CGAL_Qt4_FOUND AND QT4_FOUND )

include(${QT_USE_FILE})

# UI files (Qt Designer files)
qt4_wrap_ui( CDT_UI_FILES Polyline_simplification_2.ui )

# qrc files (resources files, that contain icons, at least)
qt4_add_resources ( CDT_RESOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/Polyline_simplification_2.qrc )

# use the Qt MOC preprocessor on classes that derives from QObject
qt4_generate_moc( ${CMAKE_CURRENT_SOURCE_DIR}/Polyline_simplification_2.cpp Polyline_simplification_2.moc )

# The executable itself.
add_executable  ( Polyline_simplification_2 ${CMAKE_CURRENT_SOURCE_DIR}/Polyline_simplification_2.cpp Polyline_simplification_2.moc ${CDT_UI_FILES} ${CDT_RESOURCE_FILES} )

# Link with Qt libraries
target_link_libraries( Polyline_simplification_2 ${QT_LIBRARIES} )

# Link with CGAL
target_link_libraries( Polyline_simplification_2 ${CGAL_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES})

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt4, and will not be compiled.")

endif()
