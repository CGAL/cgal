# This is the CMake script for compiling the AABB tree demo.

project( AABB_demo )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

# Include this package's headers first
include_directories( BEFORE ./ ./include ../../include )


# Find Qt itself
set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )


find_package(CGAL REQUIRED COMPONENTS Qt4or5)

include(${CGAL_USE_FILE})
	
find_package(Qt4or5)


# Find OpenGL
find_package(OpenGL)
find_package(QGLViewer )

if(QT${QT_VERSION_USED}_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )

  qt_wrap_ui( UI_FILES MainWindow.ui )

  include(AddFileDependencies)

  qt_generate_moc( "MainWindow.h" "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp" )
  add_file_dependencies( MainWindow_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" )

  qt_generate_moc( "Viewer.h" "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp" )
  add_file_dependencies( Viewer_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h" )

  qt_generate_moc( "Scene.h" "${CMAKE_CURRENT_BINARY_DIR}/Scene_moc.cpp" )
  add_file_dependencies( Scene_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Scene.h" )

  qt_add_resources ( RESOURCE_FILES AABB_demo.qrc )

  add_file_dependencies( AABB_demo.cpp "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp"
                                       "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp"
                                       "${CMAKE_CURRENT_BINARY_DIR}/Scene_moc.cpp" )

  add_executable  ( AABB_demo AABB_demo.cpp ${UI_FILES} ${RESOURCE_FILES} )

  # Link with Qt libraries
  target_link_libraries( AABB_demo ${QT_LIBRARIES} )

  # Link with CGAL
  target_link_libraries( AABB_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  # Link with libQGLViewer, OpenGL
  target_link_libraries( AABB_demo ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS AABB_demo )


else ( QT${QT_VERSION_USED}_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND )

  set(AABB_MISSING_DEPS "")

  if( NOT CGAL_${QT_VERSION}_FOUND )
    set(AABB_MISSING_DEPS "the CGAL ${QT_VERSION} library, ${AABB_MISSING_DEPS}")
  endif()

  if( NOT QT${QT_VERSION_USED}_FOUND )
    set(AABB_MISSING_DEPS "${QT_VERSION}, ${AABB_MISSING_DEPS}")
  endif()

  if(NOT OPENGL_FOUND)
    set(AABB_MISSING_DEPS "OpenGL, ${AABB_MISSING_DEPS}")
  endif()

  if(NOT QGLVIEWER_FOUND)
    set(AABB_MISSING_DEPS "QGLViewer, ${AABB_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${AABB_MISSING_DEPS}and will not be compiled.")

endif ( QT${QT_VERSION_USED}_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)
