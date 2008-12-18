# This is the CMake script for compiling the CGAL Polyhedron demo.

project( Polyhedron )

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

#option(POLYHEDRON_DEMO_ENABLE_FORWARD_DECL "In the Polyhedron demo, enable " OFF)
#mark_as_advanced(POLYHEDRON_DEMO_ENABLE_FORWARD_DECL)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

foreach(INCDIR ./include ../../../AABB_tree/include ../../../STL_Extension/include ../../../GraphicsView/include ../../../Surface_mesher/include ../../../filtered_kernel/include )
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${INCDIR}")
    include_directories (BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/${INCDIR}")
  endif()
endforeach()

include_directories( ./ )

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

# Find TAUCS (optionnal)
find_package(TAUCS)

if(CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )

  if(TAUCS_FOUND)
    include( ${TAUCS_USE_FILE} )
    add_definitions(-DCGAL_TAUCS_ENABLED)
  else(TAUCS_FOUND)
    message(STATUS "warning: TAUCS is not found. parametrization will not be available.")
  endif(TAUCS_FOUND)

  qt4_wrap_ui( UI_FILES MainWindow.ui )

  include(AddFileDependencies)

  qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" MainWindow_moc.cpp )
  add_file_dependencies( MainWindow_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" )

  qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h" Viewer_moc.cpp )
  add_file_dependencies( Viewer_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h" )

  qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/Scene.h" Scene_moc.cpp )
  add_file_dependencies( Scene_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Scene.h" )

  qt4_add_resources ( RESOURCE_FILES Polyhedron_3.qrc )

  if(POLYHEDRON_DEMO_ENABLE_FORWARD_DECL)
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
  else(POLYHEDRON_DEMO_ENABLE_FORWARD_DECL)
    add_file_dependencies( Polyhedron_3.cpp "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp"
                                            "${CMAKE_CURRENT_BINARY_DIR}/Scene_moc.cpp"
					    "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp" )
    add_executable  ( Polyhedron_3 Polyhedron_3.cpp ${UI_FILES} ${RESOURCE_FILES} )
  endif(POLYHEDRON_DEMO_ENABLE_FORWARD_DECL)


  # Link with Qt libraries
  target_link_libraries( Polyhedron_3 ${QT_LIBRARIES} )

  # Link with CGAL
  target_link_libraries( Polyhedron_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  # Link with libQGLViewer, OpenGL
  target_link_libraries( Polyhedron_3 ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Polyhedron_3 )

  if(TAUCS_FOUND)
    # Link with TAUCS
    target_link_libraries(Polyhedron_3 ${TAUCS_LIBRARIES})
  endif(TAUCS_FOUND)

else (CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  set(POLYHEDRON_MISSING_DEPS "")

  if(NOT CGAL_Qt4_FOUND)
    set(POLYHEDRON_MISSING_DEPS "the CGAL Qt4 library, ${POLYHEDRON_MISSING_DEPS}")
  endif()

  if(NOT QT4_FOUND)
    set(POLYHEDRON_MISSING_DEPS "Qt4, ${POLYHEDRON_MISSING_DEPS}")
  endif()

  if(NOT OPENGL_FOUND)
    set(POLYHEDRON_MISSING_DEPS "OpenGL, ${POLYHEDRON_MISSING_DEPS}")
  endif()

  if(NOT QGLVIEWER_FOUND)
    set(POLYHEDRON_MISSING_DEPS "QGLViewer, ${POLYHEDRON_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${POLYHEDRON_MISSING_DEPS}and will not be compiled.")

endif (CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)
