# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (simpleViewer)

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
find_package(Qt4)

find_package(OpenGL)
find_package(QGLViewer)

if ( CGAL_FOUND AND CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND )

  include(${QT_USE_FILE})

  include_directories (${QGLVIEWER_INCLUDE_DIR})
  include_directories (BEFORE ../../include ./ )

  # ui file, created wih Qt Designer
  qt4_wrap_ui( uis Mainwindow.ui )

  # qrc files (resources files, that contain icons, at least)
  #qt4_add_resources ( RESOURCE_FILES ./Alpha_shape_3.qrc )

  qt4_generate_moc( main.cpp Mainwindow.moc)

  add_executable  ( simpleViewer Mainwindow.moc main.cpp ${uis}  ${RESOURCE_FILES} )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS simpleViewer )

  target_link_libraries( simpleViewer ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
  target_link_libraries( simpleViewer ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES} )
  target_link_libraries( simpleViewer ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

else()

  message(STATUS "NOTICE: This demo requires CGAL, the QGLViewer, OpenGL and Qt4, and will not be compiled.")

endif()
