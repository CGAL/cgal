# This is the CMake script for compiling the CGAL Mesh_3 demo.

project( Mesh_3 )

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

#option(MESH_3_DEMO_ENABLE_FORWARD_DECL "In the Mesh_3 demo, enable " OFF)
#mark_as_advanced(MESH_3_DEMO_ENABLE_FORWARD_DECL)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# Use packages improved since CGAL 3.4
foreach(DEP_PKG AABB_tree STL_Extension GraphicsView Surface_mesher Filtered_kernel Profiling_tools Mesh_3 Triangulation_3 )
  foreach(CGAL_SVN_TRUNK ../../../ ../../trunk ..)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${CGAL_SVN_TRUNK}/${DEP_PKG}/include")
      include_directories (BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/${CGAL_SVN_TRUNK}/${DEP_PKG}/include")
      break()
    endif()
  endforeach()
endforeach()

include_directories(../../../Polyhedron/demo/Polyhedron)

# Include this package's headers first
include_directories( BEFORE ./ ./include ../../include )

# Find CGAL and CGAL Qt4
find_package(CGAL COMPONENTS Qt4)
include( ${CGAL_USE_FILE} )

# Find Qt4 itself
set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )
find_package(Qt4)

# Find OpenGL
find_package(OpenGL)

# Find QGLViewer
if(QT4_FOUND)
  include(${QT_USE_FILE})
  find_package(QGLViewer )
endif(QT4_FOUND)

# Find LAPACK (optional), for curvatures estimation
find_package(LAPACK)

if(CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )

  qt4_wrap_ui( UI_FILES MainWindow.ui )

  include(AddFileDependencies)

  qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" MainWindow_moc.cpp )
  add_file_dependencies( MainWindow_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" )

  qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h" Viewer_moc.cpp )
  add_file_dependencies( Viewer_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h" )

  qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/Scene.h" Scene_moc.cpp )
  add_file_dependencies( Scene_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Scene.h" )

  qt4_add_resources ( RESOURCE_FILES Mesh_3.qrc )

  add_definitions(-DUSE_FORWARD_DECL)
  add_definitions(-DQT_STATICPLUGIN)
    add_executable  ( Mesh_3 MainWindow.cpp
      Mesh_3.cpp
      Viewer.cpp
      Scene.cpp
      MainWindow_moc.cpp
      Scene_moc.cpp
      Viewer_moc.cpp
      ${UI_FILES} ${RESOURCE_FILES} )
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS Mesh_3 )

  # Link with Qt libraries
  target_link_libraries( Mesh_3 ${QT_LIBRARIES} )

  # Link with CGAL
  target_link_libraries( Mesh_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  # Link with libQGLViewer, OpenGL
  target_link_libraries( Mesh_3 ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

  # Link with the scene_item library.
#  target_link_libraries( Mesh_3 scene_item )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Mesh_3 )

else (CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  set(MESH_3_MISSING_DEPS "")

  if(NOT CGAL_Qt4_FOUND)
    set(MESH_3_MISSING_DEPS "the CGAL Qt4 library, ${MESH_3_MISSING_DEPS}")
  endif()

  if(NOT QT4_FOUND)
    set(MESH_3_MISSING_DEPS "Qt4, ${MESH_3_MISSING_DEPS}")
  endif()

  if(NOT OPENGL_FOUND)
    set(MESH_3_MISSING_DEPS "OpenGL, ${MESH_3_MISSING_DEPS}")
  endif()

  if(NOT QGLVIEWER_FOUND)
    set(MESH_3_MISSING_DEPS "QGLViewer, ${MESH_3_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${MESH_3_MISSING_DEPS}and will not be compiled.")

endif (CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)
