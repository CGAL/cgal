# This is the CMake script for compiling the PCA demo.

project( PCA_demo )

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

foreach(INCDIR ../../include ../../../STL_Extension/include ../../../GraphicsView/include ../../../filtered_kernel/include )
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${INCDIR}")
    include_directories (BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/${INCDIR}")
  endif()
endforeach()

include_directories( ./ )

# Find CGAL and CGAL Qt5
find_package(CGAL COMPONENTS Qt5)
include( ${CGAL_USE_FILE} )

# Find Qt5 itself
find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL)

# Find OpenGL
find_package(OpenGL)

# Find QGLViewer
if(Qt5_FOUND)
  find_package(QGLViewer )
endif(Qt5_FOUND)

if(CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )

  qt5_wrap_ui( UI_FILES MainWindow.ui )

  include(AddFileDependencies)

  qt5_generate_moc( "MainWindow.h" "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp" )
  add_file_dependencies( MainWindow_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" )

  qt5_generate_moc( "Viewer.h" "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp" )
  add_file_dependencies( Viewer_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h" )

  qt5_add_resources ( RESOURCE_FILES PCA_demo.qrc )

    add_file_dependencies( PCA_demo.cpp "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp"
                                         "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp" )
    add_executable  ( PCA_demo PCA_demo.cpp ${UI_FILES} ${RESOURCE_FILES} )
    qt5_use_modules( PCA_demo Xml Script OpenGL)
  # Link with Qt libraries
  target_link_libraries( PCA_demo ${QT_LIBRARIES} )

  # Link with CGAL
  target_link_libraries( PCA_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  # Link with libQGLViewer, OpenGL
  target_link_libraries( PCA_demo ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS PCA_demo )


else (CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  set(PCA_MISSING_DEPS "")

  if(NOT CGAL_Qt5_FOUND)
    set(PCA_MISSING_DEPS "the CGAL Qt5 library, ${PCA_MISSING_DEPS}")
  endif()

  if(NOT Qt5_FOUND)
    set(PCA_MISSING_DEPS "Qt5, ${PCA_MISSING_DEPS}")
  endif()

  if(NOT OPENGL_FOUND)
    set(PCA_MISSING_DEPS "OpenGL, ${PCA_MISSING_DEPS}")
  endif()

  if(NOT QGLVIEWER_FOUND)
    set(PCA_MISSING_DEPS "QGLViewer, ${PCA_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${PCA_MISSING_DEPS} and will not be compiled.")

endif (CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)
