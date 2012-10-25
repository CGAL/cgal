project (Circular_kernel_3)

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
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

  add_executable  ( Circular_kernel_3 Circular_kernel_3.cpp Viewer.cpp )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Circular_kernel_3 )

  target_link_libraries( Circular_kernel_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
  target_link_libraries( Circular_kernel_3 ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES} )
  target_link_libraries( Circular_kernel_3 ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

else()

  message(STATUS "NOTICE: This demo requires CGAL, the QGLViewer, OpenGL and Qt4, and will not be compiled.")

endif()
