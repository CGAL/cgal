# This is the CMake script for compiling the CGAL Point Set demo.

project( Point_set_demo )

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

#option(POINT_SET_DEMO_ENABLE_FORWARD_DECL "In the Point Set demo, enable " OFF)
#mark_as_advanced(POINT_SET_DEMO_ENABLE_FORWARD_DECL)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

foreach(INCDIR ../../../../AABB_tree/include ../../../../STL_Extension/include ../../../../GraphicsView/include ../../../../Surface_mesher/include ../../../../Filtered_kernel/include ../../../../Profiling_tools/include)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${INCDIR}")
    include_directories (BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/${INCDIR}")
  endif()
endforeach()

include_directories (BEFORE ../../../../Point_set_processing_3/include/)
include_directories( BEFORE ./ ./include ../../../include )

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

# Find LAPACK (optional), for Jet Fitting
find_package(LAPACK)

# Find TAUCS (optionnal), for Poisson
find_package(TAUCS)

if(CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )

  if(LAPACK_FOUND)
    include( ${LAPACK_USE_FILE} )
    add_definitions(-DCGAL_LAPACK_ENABLED)
  else(LAPACK_FOUND)
    message(STATUS "warning: LAPACK is not found. Jet Fitting will not be available.")
  endif(LAPACK_FOUND)

  if(TAUCS_FOUND)
    include( ${TAUCS_USE_FILE} )
    add_definitions(-DCGAL_TAUCS_ENABLED)
  else(TAUCS_FOUND)
    message(STATUS "warning: TAUCS is not found. Poisson reconstruction will not be available.")
  endif(TAUCS_FOUND)

  qt4_wrap_ui( UI_FILES MainWindow.ui )

  include(AddFileDependencies)

  qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" MainWindow_moc.cpp )
  add_file_dependencies( MainWindow_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" )

  qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h" Viewer_moc.cpp )
  add_file_dependencies( Viewer_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h" )

  qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/Scene.h" Scene_moc.cpp )
  add_file_dependencies( Scene_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Scene.h" )

  qt4_add_resources ( RESOURCE_FILES Point_set_demo.qrc )

  qt4_automoc(Scene_item.cpp
              Scene_plane_item.cpp
              Point_set_scene_item.cpp 
              Scene_polyhedron_item.cpp)

  # AUXILIARY LIBRARIES
  add_library(scene_item SHARED
    Scene_item.cpp Scene_item.moc
    Scene_item_with_display_list.cpp
    Polyhedron_demo_plugin_helper.cpp)

  add_library(scene_basic_objects SHARED
    Scene_plane_item.cpp Scene_plane_item.moc)
  target_link_libraries(scene_basic_objects scene_item ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})

  add_library(scene_polyhedron_item SHARED
    Scene_polyhedron_item.cpp Scene_polyhedron_item.moc)
  target_link_libraries(scene_polyhedron_item scene_item)

  add_library(point_set SHARED
    Point_set_scene_item.cpp Point_set_scene_item.moc)
  target_link_libraries(point_set scene_item)

  add_definitions(-DQT_STATICPLUGIN)
#  if(POINT_SET_DEMO_ENABLE_FORWARD_DECL)
    add_definitions(-DUSE_FORWARD_DECL)
    add_executable  ( Point_set_demo MainWindow.cpp
      Point_set_demo.cpp
      Viewer.cpp
      Scene.cpp
      MainWindow_moc.cpp
      Scene_moc.cpp
      Viewer_moc.cpp
      ${UI_FILES} ${RESOURCE_FILES} )
    if(TAUCS_FOUND)
#      add_executable( Point_set_demo Scene_tex_rendering.cpp Scene_tex_polyhedron_operations.cpp )
    endif()
#  else(POINT_SET_DEMO_ENABLE_FORWARD_DECL)
#    add_file_dependencies( Point_set_demo.cpp "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp"
#                                            "${CMAKE_CURRENT_BINARY_DIR}/Scene_moc.cpp"
#					    "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp" )
#    add_executable  ( Point_set_demo Point_set_demo.cpp ${UI_FILES} ${RESOURCE_FILES} )
#  endif(POINT_SET_DEMO_ENABLE_FORWARD_DECL)


  # Link with Qt libraries
  target_link_libraries( Point_set_demo ${QT_LIBRARIES} )

  # Link with CGAL
  target_link_libraries( Point_set_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  # Link with libQGLViewer, OpenGL
  target_link_libraries( Point_set_demo ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

  # Link with the scene_item library.
#  target_link_libraries( Point_set_demo scene_item )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Point_set_demo )


  ###########
  # PLUGINS #
  ###########
  remove_definitions(-DQT_STATICPLUGIN)

  macro(polyhedron_demo_plugin plugin_name plugin_implementation_base_name)
    list_split(option ARGN_TAIL ${ARGN} )
    if(NOT ${option} STREQUAL "EXCLUDE_FROM_ALL")
      set(other_sources ${ARGN})
      set(option "")
    else()
      set(other_sources ${ARGN_TAIL})
    endif()
    qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/${plugin_implementation_base_name}.cpp" ${plugin_implementation_base_name}.moc )
    add_file_dependencies( ${plugin_implementation_base_name}.moc "${CMAKE_CURRENT_SOURCE_DIR}/${plugin_implementation_base_name}.cpp" )

    add_library(${plugin_name} MODULE ${option} ${plugin_implementation_base_name}.moc ${plugin_implementation_base_name}.cpp ${other_sources})
    # Link with Qt
    target_link_libraries( ${plugin_name} ${QT_LIBRARIES} )
    # Link with scene_item
    target_link_libraries( ${plugin_name} scene_item)
    # Link with CGAL
    target_link_libraries( ${plugin_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
  endmacro(polyhedron_demo_plugin)

  polyhedron_demo_plugin(remeshing_plugin Point_set_demo_remeshing_plugin Point_set_demo_remeshing_plugin_cgal_code.cpp)
  target_link_libraries(remeshing_plugin scene_polyhedron_item)

  polyhedron_demo_plugin(inside_out_plugin Polyhedron_demo_inside_out_plugin)
  target_link_libraries(inside_out_plugin scene_polyhedron_item point_set)

  polyhedron_demo_plugin(off_plugin Point_set_demo_off_plugin)
  target_link_libraries(off_plugin scene_polyhedron_item point_set)

  polyhedron_demo_plugin(xyz_plugin Point_set_demo_xyz_plugin)
  target_link_libraries(xyz_plugin scene_polyhedron_item point_set)

  polyhedron_demo_plugin(convex_hull_plugin Point_set_demo_convex_hull_plugin)
  target_link_libraries(convex_hull_plugin scene_polyhedron_item)

  polyhedron_demo_plugin(pca_plugin Point_set_demo_pca_plugin)
  target_link_libraries(pca_plugin scene_polyhedron_item scene_basic_objects)

  polyhedron_demo_plugin(mesh_simplification_plugin Point_set_demo_mesh_simplification_plugin)
  target_link_libraries(mesh_simplification_plugin scene_polyhedron_item)

  polyhedron_demo_plugin(cut_plugin Polyhedron_demo_cut_plugin)
  target_link_libraries(cut_plugin scene_polyhedron_item scene_basic_objects)

else (CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  set(POINT_SET_DEMO_MISSING_DEPS "")

  if(NOT CGAL_Qt4_FOUND)
    set(POINT_SET_DEMO_MISSING_DEPS "the CGAL Qt4 library, ${POINT_SET_DEMO_MISSING_DEPS}")
  endif()

  if(NOT QT4_FOUND)
    set(POINT_SET_DEMO_MISSING_DEPS "Qt4, ${POINT_SET_DEMO_MISSING_DEPS}")
  endif()

  if(NOT OPENGL_FOUND)
    set(POINT_SET_DEMO_MISSING_DEPS "OpenGL, ${POINT_SET_DEMO_MISSING_DEPS}")
  endif()

  if(NOT QGLVIEWER_FOUND)
    set(POINT_SET_DEMO_MISSING_DEPS "QGLViewer, ${POINT_SET_DEMO_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${POINT_SET_DEMO_MISSING_DEPS}and will not be compiled.")

endif (CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)
