# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Triangulation_demo_3)

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

#New for Qt5 version !
option (USE_QT5 "Use Qt5 version instead of Qt4" ON)

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )

if(USE_QT5)
	find_package(CGAL COMPONENTS Qt5)
	include(${CGAL_USE_FILE})

	find_package(Qt5)

else(USE_QT5)
	find_package(CGAL  COMPONENTS Qt4)
	include(${CGAL_USE_FILE})

	find_package(Qt4)
	if(QT4_FOUND)
		message("Qt4 found")
		set(QT4 TRUE) 
                include( ${QT_USE_FILE} )
	endif(QT4_FOUND)
endif (USE_QT5)

include(QtChoice)

find_package(OpenGL)
find_package(QGLViewer)

if ( CGAL_FOUND AND (QT4 OR QT5) AND OPENGL_FOUND AND QGLVIEWER_FOUND )

  include_directories (${QGLVIEWER_INCLUDE_DIR})
  include_directories (BEFORE ../../include ./ )

  FILE (GLOB CPP *.cpp)

  # ui files, created with Qt Designer
  qt_wrap_ui( uis MainWindow.ui )

  # qrc files (resources files, that contain icons, at least)
  qt_add_resources ( RESOURCE_FILES ./T3_demo.qrc )

  # moc files (meta object compiler files, run moc preprocessor on the files with Widget defs)
  qt_automoc( MainWindow.cpp Viewer.cpp PreferenceDlg.cpp )


  # cpp files
  add_executable  ( T3_demo
			${CPP} 
			${uis}  
			${RESOURCE_FILES} )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS T3_demo )

  target_link_libraries( T3_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
  target_link_libraries( T3_demo ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES} )
  target_link_libraries( T3_demo ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

else( CGAL_FOUND AND (QT4 OR QT5) AND OPENGL_FOUND AND QGLVIEWER_FOUND )

  set(TRIANGULATION_3_MISSING_DEPS "")

  if(NOT CGAL_FOUND)
    set(TRIANGULATION_3_MISSING_DEPS "the CGAL Qt library, ${TRIANGULATION_3_MISSING_DEPS}")
  endif()      	   

  if(NOT QT4)
    set(TRIANGULATION_3_MISSING_DEPS "Qt4, ${TRIANGULATION_3_MISSING_DEPS}")
  endif() 

  if(NOT QT5)
    set(TRIANGULATION_3_MISSING_DEPS "Qt5, ${TRIANGULATION_3_MISSING_DEPS}")  
  endif()     	   

  if(NOT OPENGL_FOUND)
    set(TRIANGULATION_3_MISSING_DEPS "OpenGL, ${TRIANGULATION_3_MISSING_DEPS}")
  endif()      	   

  if(NOT QGLVIEWER_FOUND)
    set(TRIANGULATION_3_MISSING_DEPS "QGLViewer, ${TRIANGULATION_3_MISSING_DEPS}")
  endif()      	   


  message(STATUS "NOTICE: This demo requires ${TRIANGULATION_3_MISSING_DEPS}and will not be compiled.")

endif( CGAL_FOUND AND (QT4 OR QT5) AND OPENGL_FOUND AND QGLVIEWER_FOUND )
