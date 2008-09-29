# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Alpha_shape_3)

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL REQUIRED  Qt4)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALQt4Program)

  set( QT_USE_QTOPENGL TRUE )
  set( QT_USE_QTXML    TRUE )
  set( QT_USE_QTMAIN   TRUE )
  set( QT_USE_QTSCRIPT  TRUE )

  find_package(QT4       REQUIRED)
  include(${QT_USE_FILE})
  find_package(OpenGL    REQUIRED)
  find_package(QGLViewer REQUIRED)

  include_directories (${QGLVIEWER_INCLUDE_DIR})
 include_directories (BEFORE ../../include)

  qt4_wrap_ui( uis MainWindow.ui )

  qt4_automoc( MainWindow.cpp Viewer.cpp)

  add_executable  ( Alpha_shape_3 Alpha_shape_3.cpp MainWindow.cpp Viewer.cpp ${uis} )
  target_link_libraries( Alpha_shape_3 ${QT_LIBRARIES} ${QGLVIEWER_LIBRARY} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )
endif()
