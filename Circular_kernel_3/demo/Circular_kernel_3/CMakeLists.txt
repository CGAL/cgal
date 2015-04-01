project (Circular_kernel_3)

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

find_package(CGAL REQUIRED COMPONENTS Qt5)
include(${CGAL_USE_FILE})

find_package(Qt5 COMPONENTS Xml Script OpenGL)

find_package(OpenGL)
find_package(QGLViewer)

if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND )


  include_directories (${QGLVIEWER_INCLUDE_DIR})
  include_directories (BEFORE ../../include ./ )

  add_executable  ( Circular_kernel_3 Circular_kernel_3.cpp Viewer.cpp )

  qt5_use_modules( Circular_kernel_3 Xml Script OpenGL)

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Circular_kernel_3 )

  target_link_libraries( Circular_kernel_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
  target_link_libraries( Circular_kernel_3 ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES} )
  target_link_libraries( Circular_kernel_3 ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

else()

  message(STATUS "NOTICE: This demo requires CGAL, the QGLViewer, OpenGL and Qt5, and will not be compiled.")

endif()
