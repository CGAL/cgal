# This is the CMake script for compiling the CGAL Mesh_3 demo.

# Add directory containing implicit function plugin source files
add_subdirectory("implicit_functions")

project( Mesh_3_demo )

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

#option(MESH_3_DEMO_ENABLE_FORWARD_DECL "In the Mesh_3 demo, enable " OFF)
#mark_as_advanced(MESH_3_DEMO_ENABLE_FORWARD_DECL)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# Use packages improved since CGAL 3.4
#foreach(DEP_PKG AABB_tree STL_Extension GraphicsView Surface_mesher Filtered_kernel Profiling_tools Mesh_3 Triangulation_3 )
#  foreach(CGAL_SVN_TRUNK ../../../ ../../trunk ..)
#    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${CGAL_SVN_TRUNK}/${DEP_PKG}/include")
#      include_directories (BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/${CGAL_SVN_TRUNK}/${DEP_PKG}/include")
#      break()
#    endif()
#  endforeach()
#endforeach()

#include_directories(../../../Polyhedron/demo/Polyhedron)
#link_directories(../../../Polyhedron/demo/Polyhedron)

# Include this package's headers first
include_directories( BEFORE ./ ./include ../../include )

#include_directories(BEFORE ../../../../candidate-packages/Mesh_3/include/)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR} )

# Find CGAL and CGAL Qt4
find_package(CGAL COMPONENTS Qt4 ImageIO)
include( ${CGAL_USE_FILE} )

# Find Qt4 itself
set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )
find_package(Qt4)

# Find OpenGL
find_package(OpenGL)

# Find GLEW (for OpenGL-1.5 and OpenGL extensions)
find_package(GLEW)


# Find QGLViewer
if(QT4_FOUND)
  include(${QT_USE_FILE})
  find_package(QGLViewer )
endif(QT4_FOUND)


if(CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  if(GLEW_FOUND)
    include_directories ( ${GLEW_INCLUDE_DIR} )
    add_definitions(-DSCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE)
  else(GLEW_FOUND)
    message(STATUS "NOTICE: GLEW library is not found. 3D images rendering will not be available.")
  endif(GLEW_FOUND)


  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )

  qt4_wrap_ui( UI_FILES ui_files/MainWindow.ui )
  qt4_wrap_ui( meshingUI_FILES  ui_files/Meshing_dialog.ui)
  qt4_wrap_ui( optimUI_FILES  ui_files/Smoother_dialog.ui ui_files/LocalOptim_dialog.ui)
  qt4_wrap_ui( ribUI_FILES  ui_files/Rib_dialog.ui)
  qt4_wrap_ui( funcUI_FILES  ui_files/Function_dialog.ui)

  include(AddFileDependencies)

  # Shared demo files directory
  set(DEMO_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/CGAL_demo")
  set(DEMO_H_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include/CGAL_demo")

  qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" MainWindow_moc.cpp )
  add_file_dependencies( MainWindow_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" )

  qt4_generate_moc( "${DEMO_H_DIR}/Viewer.h" Viewer_moc.cpp )
  add_file_dependencies( Viewer_moc.cpp "${DEMO_H_DIR}/Viewer.h" )

  qt4_generate_moc( "${DEMO_H_DIR}/Scene.h" Scene_moc.cpp )
  add_file_dependencies( Scene_moc.cpp "${DEMO_H_DIR}/Scene.h" )

  qt4_generate_moc( "${DEMO_H_DIR}/Scene_item.h" Scene_item_moc.cpp )
  add_file_dependencies( Scene_item_moc.cpp "${DEMO_H_DIR}/Scene_item.h" )

  qt4_add_resources ( RESOURCE_FILES Mesh_3.qrc )

  qt4_automoc( Scene_polyhedron_item.cpp )
  qt4_automoc( Scene_polygon_soup.cpp )
  qt4_automoc( Scene_segmented_image_item.cpp )
  qt4_automoc( Scene_c3t3_item.cpp )
  qt4_automoc( Scene_implicit_function_item.cpp )
  qt4_automoc( Optimizer_thread.cpp )
  qt4_automoc( Meshing_thread.cpp )

  # put plugins (which are shared libraries) at the same location as
  # executable files
  set(LIBRARY_OUTPUT_PATH ${RUNTIME_OUTPUT_PATH})

  # AUXILIARY LIBRARIES
  add_library(viewer SHARED
    ${DEMO_SRC_DIR}/Viewer.cpp Viewer_moc.cpp)
  target_link_libraries(viewer ${QGLVIEWER_LIBRARIES})

  add_library(scene_item SHARED
    ${DEMO_SRC_DIR}/Scene_item.cpp
    Scene_item_moc.cpp
    ${DEMO_SRC_DIR}/Scene_item_with_display_list.cpp
    ${DEMO_SRC_DIR}/Plugin_helper.cpp)

  add_library(scene_segmented_image_item SHARED
    Scene_segmented_image_item.cpp Scene_segmented_image_item.moc)
  target_link_libraries(scene_segmented_image_item scene_item)

  if(GLEW_FOUND)
    target_link_libraries(scene_segmented_image_item ${GLEW_LIBRARIES})
  endif()

  add_library(scene_polyhedron_item SHARED
    Scene_polyhedron_item.cpp Scene_polyhedron_item.moc)
  target_link_libraries(scene_polyhedron_item scene_item)

  add_library(polygon_soup SHARED
    Scene_polygon_soup.cpp Scene_polygon_soup.moc)
  target_link_libraries(polygon_soup scene_item)

  add_library(scene_c3t3_item SHARED
    Scene_c3t3_item.cpp Scene_c3t3_item.moc)
  target_link_libraries(scene_c3t3_item scene_item ${QGLVIEWER_LIBRARIES} )

  add_library(scene_implicit_function_item SHARED
    Scene_implicit_function_item.cpp Scene_implicit_function_item.moc Color_ramp.cpp)
  target_link_libraries(scene_implicit_function_item scene_item ${QGLVIEWER_LIBRARIES} )


  add_definitions(-DUSE_FORWARD_DECL)
  add_definitions(-DQT_STATICPLUGIN)
    add_executable  ( Mesh_3 MainWindow.cpp
      Mesh_3.cpp
      ${DEMO_SRC_DIR}/Scene.cpp
      MainWindow_moc.cpp
      Scene_moc.cpp
      ${UI_FILES} ${RESOURCE_FILES} )
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS Mesh_3 )

  # Link with Qt libraries
  target_link_libraries( Mesh_3 ${QT_LIBRARIES} )

  # Link with CGAL
  target_link_libraries( Mesh_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  # Link with libQGLViewer, OpenGL
  target_link_libraries( Mesh_3 ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

#  # Link with the scene_item library.
#  target_link_libraries( Mesh_3 scene_item )

  target_link_libraries( Mesh_3 viewer )


  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Mesh_3 )


  ###########
  # PLUGINS #
  ###########
  remove_definitions(-DQT_STATICPLUGIN)

  # polyhedron_demo_plugin is a macro copied from demo/Polyhedron/CMakeLists.txt
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
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${plugin_name} )
    # Link with Qt
    target_link_libraries( ${plugin_name} ${QT_LIBRARIES} )
    # Link with scene_item
    target_link_libraries( ${plugin_name} scene_item)
    # Link with CGAL
    target_link_libraries( ${plugin_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
  endmacro(polyhedron_demo_plugin)

  polyhedron_demo_plugin(io_image_plugin Io_image_plugin)
  target_link_libraries(io_image_plugin scene_segmented_image_item)

  if(GLEW_FOUND)
    target_link_libraries(io_image_plugin ${GLEW_LIBRARIES})
  endif()

  polyhedron_demo_plugin(io_off_plugin Io_off_plugin)
  target_link_libraries(io_off_plugin scene_polyhedron_item polygon_soup)

  polyhedron_demo_plugin(io_c3t3_plugin Io_c3t3_plugin)
  target_link_libraries(io_c3t3_plugin scene_c3t3_item)

  if ( EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/c3t3_rib_exporter_plugin.cpp" )
    polyhedron_demo_plugin(c3t3_rib_exporter_plugin C3t3_rib_exporter_plugin ${ribUI_FILES})
    target_link_libraries(c3t3_rib_exporter_plugin scene_c3t3_item viewer)
  endif()

  polyhedron_demo_plugin(io_implicit_function_plugin Io_implicit_function_plugin ${funcUI_FILES})
  target_link_libraries(io_implicit_function_plugin scene_implicit_function_item)

  polyhedron_demo_plugin(mesh_3_plugin Mesh_3_plugin 
                         Mesh_3_plugin_polyhedron_cgal_code.cpp
                         Mesh_3_plugin_image_cgal_code.cpp
                         Mesh_3_plugin_implicit_function_cgal_code.cpp
                         Meshing_thread.cpp
                         Scene_c3t3_item.moc
                         ${meshingUI_FILES})

  target_link_libraries(mesh_3_plugin
                        scene_c3t3_item
                        scene_polyhedron_item
                        scene_segmented_image_item
                        scene_implicit_function_item
                        ${QGLVIEWER_LIBRARIES}
                        ${OPENGL_gl_LIBRARY}
                        ${OPENGL_glu_LIBRARY})

  polyhedron_demo_plugin(mesh_3_optimization_plugin Mesh_3_optimization_plugin 
                         Mesh_3_optimization_plugin_cgal_code.cpp
                         Optimizer_thread.cpp
                         Scene_c3t3_item.moc
                         ${optimUI_FILES})

  target_link_libraries(mesh_3_optimization_plugin
                        scene_c3t3_item
                        scene_polyhedron_item
                        scene_segmented_image_item
                        scene_implicit_function_item
                        ${QGLVIEWER_LIBRARIES}
                        ${OPENGL_gl_LIBRARY}
                        ${OPENGL_glu_LIBRARY})


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


