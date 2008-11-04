# This is the CMake script for compiling a CGAL application.

project( Polyhedron )

cmake_minimum_required(VERSION 2.4)

option(ENABLE_FORWARD_DECL "Enable the forward declaration of Polyhedron and Nef_polyhedron." OFF)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

foreach(INCDIR ./include ../../../AABB_tree/include ../../../STL_Extension/include ../../../GraphicsView/include ../../../Surface_mesher/include ../../../filtered_kernel/include )
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${INCDIR})
    include_directories (BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/${INCDIR})
  endif()
endforeach()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../../../GraphicsView/src/CGALQt4)
  link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../GraphicsView/src/CGALQt4)
endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../../../GraphicsView/src/CGALQt4)

# Find CGAL and CGAL Qt4
find_package(CGAL QUIET COMPONENTS Qt4)
include( ${CGAL_USE_FILE} )

# Find Boost and Boost.Thread
find_package(Boost 1.33.1 QUIET COMPONENTS thread)

# Find Qt4 itself
set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )
find_package(Qt4 QUIET)

# Find OpenGL
find_package(OpenGL QUIET)

# Find QGLViewer
if(QT_FOUND)
  include(${QT_USE_FILE})
  find_package(QGLViewer QUIET )
endif(QT_FOUND)

# Link with BLAS, LAPACK and TAUCS (optional)
include( CGAL_UseTAUCS )
if (NOT TAUCS_FOUND)
  message(STATUS "warning: TAUCS is not found. parameterization will not be available.")
endif()

if(CGAL_Qt4_FOUND AND Boost_FOUND AND Boost_THREAD_FOUND AND QT_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  qt4_wrap_ui( UI_FILES MainWindow.ui )

  include(AddFileDependencies)

  qt4_generate_moc( MainWindow.h MainWindow_moc.cpp )
  add_file_dependencies(  MainWindow_moc.cpp MainWindow.h )

  qt4_generate_moc( Viewer.h Viewer_moc.cpp )
  add_file_dependencies( Viewer_moc.cpp Viewer.h )

  qt4_generate_moc( Scene.h Scene_moc.cpp )
  add_file_dependencies( Scene_moc.cpp Scene.h )

  qt4_add_resources ( RESOURCE_FILES Polyhedron_3.qrc )

  if(ENABLE_FORWARD_DECL)
    add_definitions(-DUSE_FORWARD_DECL)
    add_executable  ( Polyhedron_3 MainWindow.cpp
      Polyhedron_3.cpp
      Viewer.cpp
      Scene.cpp
      Scene_nef_polyhedron_operations.cpp
      Scene_polyhedron_operations.cpp
      Scene_nef_and_polyhedron_operations.cpp
      Scene_rendering.cpp
      Scene_nef_rendering.cpp
      MainWindow_boolean_operations.cpp
      MainWindow_convex_hull.cpp
      MainWindow_curvature_estimation.cpp
      MainWindow_inside_out.cpp
      MainWindow_kernel.cpp
      MainWindow_pca.cpp
      MainWindow_remeshing.cpp
      MainWindow_self_intersection.cpp
      MainWindow_simplify.cpp
      MainWindow_subdivision_methods.cpp
      MainWindow_parameterization.cpp
      MainWindow_moc.cpp
      Scene_moc.cpp
      Viewer_moc.cpp
      texture.cpp
      ${UI_FILES} ${RESOURCE_FILES} )
    if(TAUCS_FOUND)
      add_executable( Polyhedron_3 Scene_tex_rendering.cpp Scene_tex_polyhedron_operations.cpp )
    endif()
  else(ENABLE_FORWARD_DECL)
    add_file_dependencies( Polyhedron_3.cpp MainWindow_moc.cpp Scene_moc.cpp Viewer_moc.cpp )
    add_executable  ( Polyhedron_3 Polyhedron_3.cpp ${UI_FILES} ${RESOURCE_FILES} )
  endif(ENABLE_FORWARD_DECL)


  # Link with Qt libraries
  target_link_libraries( Polyhedron_3 ${QT_LIBRARIES} )

  # Link with CGAL
  target_link_libraries( Polyhedron_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  # Link with libQGLViewer, OpenGL
  target_link_libraries( Polyhedron_3 ${QGLVIEWER_LIBRARY} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

else (CGAL_Qt4_FOUND AND Boost_FOUND AND Boost_THREAD_FOUND AND QT_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  set(POLYHEDRON_MISSING_DEPS "")

  if(NOT CGAL_Qt4_FOUND)
    set(POLYHEDRON_MISSING_DEPS "the CGAL Qt4 library, ${POLYHEDRON_MISSING_DEPS}")
  endif()

  if(NOT Boost_THREAD_FOUND)
    set(POLYHEDRON_MISSING_DEPS "Boost Thread, ${POLYHEDRON_MISSING_DEPS}")
  endif()

  if(NOT Boost_FOUND)
    set(POLYHEDRON_MISSING_DEPS "Boost, ${POLYHEDRON_MISSING_DEPS}")
  endif()

  if(NOT QT_FOUND)
    set(POLYHEDRON_MISSING_DEPS "Qt4, ${POLYHEDRON_MISSING_DEPS}")
  endif()

  if(NOT OPENGL_FOUND)
    set(POLYHEDRON_MISSING_DEPS "OpenGL, ${POLYHEDRON_MISSING_DEPS}")
  endif()

  if(NOT QGLVIEWER_FOUND)
    set(POLYHEDRON_MISSING_DEPS "QGLViewer, ${POLYHEDRON_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${POLYHEDRON_MISSING_DEPS}and will not be compiled.")

endif (CGAL_Qt4_FOUND AND Boost_FOUND AND Boost_THREAD_FOUND AND QT_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)
