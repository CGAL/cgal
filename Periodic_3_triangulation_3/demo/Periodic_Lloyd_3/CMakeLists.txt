# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Periodic_Lloyd_3)

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL  COMPONENTS Qt4)
include(${CGAL_USE_FILE})

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )
set( QT_USE_QTASSISTANT TRUE )
find_package(Qt4)

find_package(OpenGL)

if(QT4_FOUND)
  include( ${QT_USE_FILE} )
  find_package(QGLViewer)
endif(QT4_FOUND)

if ( CGAL_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND )

  include(${QT_USE_FILE})

  include_directories (${QGLVIEWER_INCLUDE_DIR})
  include_directories (BEFORE ../../include ./ )

  # ui file, created wih Qt Designer
  qt4_wrap_ui( uis MainWindow.ui )

  # qrc files (resources files, that contain icons, at least)
  qt4_add_resources ( RESOURCE_FILES ./Periodic_Lloyd_3.qrc )

  qt4_automoc( MainWindow.cpp Viewer.cpp)

  add_executable  ( Periodic_Lloyd_3 Periodic_Lloyd_3.cpp MainWindow.cpp Viewer.cpp ${uis}  ${RESOURCE_FILES} )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Periodic_Lloyd_3 )

  target_link_libraries( Periodic_Lloyd_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
  target_link_libraries( Periodic_Lloyd_3 ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES} )
  target_link_libraries( Periodic_Lloyd_3 ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

else( CGAL_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND )

  set(PERIODIC_LLOYD_MISSING_DEPS "")

  if(NOT CGAL_FOUND)
    set(PERIODIC_LLOYD_MISSING_DEPS "the CGAL Qt4 library, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()      	   

  if(NOT QT4_FOUND)
    set(PERIODIC_LLOYD_MISSING_DEPS "Qt4, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()      	   

  if(NOT OPENGL_FOUND)
    set(PERIODIC_LLOYD_MISSING_DEPS "OpenGL, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()      	   

  if(NOT QGLVIEWER_FOUND)
    set(PERIODIC_LLOYD_MISSING_DEPS "QGLViewer, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()      	   


  message(STATUS "NOTICE: This demo requires ${PERIODIC_LLOYD_MISSING_DEPS}and will not be compiled.")

endif( CGAL_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND )
