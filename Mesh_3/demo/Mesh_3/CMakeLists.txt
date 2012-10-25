# This is the CMake script for compiling the CGAL Mesh_3 demo.

project( Mesh_3_demo )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

# Let plugins be compiled in the same directory as the executable.
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

#option(MESH_3_DEMO_ENABLE_FORWARD_DECL "In the Mesh_3 demo, enable " OFF)
#mark_as_advanced(MESH_3_DEMO_ENABLE_FORWARD_DECL)

# Include this package's headers first
include_directories( BEFORE ./ ./include ../../include )

# Disable deprecated Mesh_3 names
add_definitions(-DCGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
                -DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)


# Add specific Find.cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules )

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

# Find everybodies darling
find_package(Boost COMPONENTS thread system)

# Find QGLViewer
if(QT4_FOUND)
  include(${QT_USE_FILE})
  find_package(QGLViewer )
endif(QT4_FOUND)


if(CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND AND Boost_FOUND)

  # Add directory containing implicit function plugin source files
  add_subdirectory("implicit_functions")

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
  qt4_wrap_ui( imgUI_FILES ui_files/Image_res_dialog.ui)

  include(AddFileDependencies)

  # Shared demo files directory
  set(DEMO_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/CGAL_demo")
  set(DEMO_H_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include/CGAL_demo")

  qt4_generate_moc( "MainWindow.h" "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp" )
  add_file_dependencies( MainWindow_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" )

  qt4_generate_moc( "${DEMO_H_DIR}/Viewer.h" Viewer_moc.cpp )
  add_file_dependencies( Viewer_moc.cpp "${DEMO_H_DIR}/Viewer.h" )

  qt4_generate_moc( "${DEMO_H_DIR}/Scene.h" Scene_moc.cpp )
  add_file_dependencies( Scene_moc.cpp "${DEMO_H_DIR}/Scene.h" )

  qt4_generate_moc( "${DEMO_H_DIR}/Scene_item.h" Scene_item_moc.cpp )
  add_file_dependencies( Scene_item_moc.cpp "${DEMO_H_DIR}/Scene_item.h" )

  qt4_wrap_cpp( VOLUME_MOC_OUTFILES ${CMAKE_CURRENT_SOURCE_DIR}/Volume_plane_thread.h )
  qt4_wrap_cpp( VOLUME_MOC_OUTFILES ${CMAKE_CURRENT_SOURCE_DIR}/Volume_plane_interface.h )

  qt4_add_resources ( RESOURCE_FILES Mesh_3.qrc )

  qt4_automoc( Scene_polyhedron_item.cpp )
  qt4_automoc( Scene_polygon_soup.cpp )
  qt4_automoc( Scene_segmented_image_item.cpp )
  qt4_automoc( Scene_c3t3_item.cpp )
  qt4_automoc( Scene_implicit_function_item.cpp )
  qt4_automoc( Optimizer_thread.cpp )
  qt4_automoc( Meshing_thread.cpp )
  qt4_automoc( Volume_plane_intersection.cpp )

  # put plugins (which are shared libraries) at the same location as
  # executable files
  set(LIBRARY_OUTPUT_PATH ${RUNTIME_OUTPUT_PATH})

  # prefix for Mesh_3 libs to avoid conflicts with Polyhedron demo
  set(MESH_3_LIB_PREFIX "mesh_3_demo_")

  # AUXILIARY LIBRARIES
  set(VIEWER_LIB "${MESH_3_LIB_PREFIX}viewer")
  add_library(${VIEWER_LIB} SHARED
    ${DEMO_SRC_DIR}/Viewer.cpp Viewer_moc.cpp)
  target_link_libraries(${VIEWER_LIB} ${QGLVIEWER_LIBRARIES} ${QT_LIBRARIES} ${OPENGL_LIBRARIES})
  set_target_properties(${VIEWER_LIB} PROPERTIES DEFINE_SYMBOL viewer_EXPORTS)

  set(SCENE_ITEM_LIB "${MESH_3_LIB_PREFIX}scene_item")
  add_library(${SCENE_ITEM_LIB} SHARED
    ${DEMO_SRC_DIR}/Scene_item.cpp
    Scene_item_moc.cpp
    ${DEMO_SRC_DIR}/Scene_item_with_display_list.cpp
    ${DEMO_SRC_DIR}/Plugin_helper.cpp)
  set_target_properties(${SCENE_ITEM_LIB} PROPERTIES DEFINE_SYMBOL scene_item_EXPORTS)
  target_link_libraries(${SCENE_ITEM_LIB} ${QT_LIBRARIES} ${OPENGL_LIBRARIES})

  set(SCENE_SEGMENTED_IMAGE_ITEM_LIB "${MESH_3_LIB_PREFIX}scene_segmented_image_item")
  add_library(${SCENE_SEGMENTED_IMAGE_ITEM_LIB} SHARED
    Scene_segmented_image_item.cpp Scene_segmented_image_item.moc)
  target_link_libraries(${SCENE_SEGMENTED_IMAGE_ITEM_LIB} ${SCENE_ITEM_LIB} ${CGAL_LIBRARIES})
  set_target_properties(${SCENE_SEGMENTED_IMAGE_ITEM_LIB} PROPERTIES DEFINE_SYMBOL scene_segmented_image_item_EXPORTS)

  if(GLEW_FOUND)
    target_link_libraries(${SCENE_SEGMENTED_IMAGE_ITEM_LIB} ${GLEW_LIBRARIES})
  endif()

  set(SCENE_POLYHEDRON_ITEM_LIB "${MESH_3_LIB_PREFIX}scene_polyhedron_item")
  add_library(${SCENE_POLYHEDRON_ITEM_LIB} SHARED
    Scene_polyhedron_item.cpp Scene_polyhedron_item.moc)
  target_link_libraries(${SCENE_POLYHEDRON_ITEM_LIB} ${SCENE_ITEM_LIB} ${CGAL_LIBRARIES})
  set_target_properties(${SCENE_POLYHEDRON_ITEM_LIB} PROPERTIES DEFINE_SYMBOL scene_polyhedron_item_EXPORTS)

  set(POLYGON_SOUP_LIB "${MESH_3_LIB_PREFIX}polygon_soup")
  add_library(${POLYGON_SOUP_LIB} SHARED
    Scene_polygon_soup.cpp Scene_polygon_soup.moc)
  target_link_libraries(${POLYGON_SOUP_LIB} ${SCENE_ITEM_LIB} ${CGAL_LIBRARIES})
  set_target_properties(${POLYGON_SOUP_LIB} PROPERTIES DEFINE_SYMBOL polygon_soup_EXPORTS)

  set(SCENE_C3T3_ITEM_LIB "${MESH_3_LIB_PREFIX}scene_c3t3_item")
  add_library(${SCENE_C3T3_ITEM_LIB} SHARED
    Scene_c3t3_item.cpp Scene_c3t3_item.moc)
  target_link_libraries(${SCENE_C3T3_ITEM_LIB} ${SCENE_ITEM_LIB} ${QGLVIEWER_LIBRARIES} ${QT_LIBRARIES} ${CGAL_LIBRARIES} ${Boost_LIBRARIES})
  set_target_properties(${SCENE_C3T3_ITEM_LIB} PROPERTIES DEFINE_SYMBOL scene_c3t3_item_EXPORTS)

  set(SCENE_IMPLICIT_FUNCTION_ITEM_LIB "${MESH_3_LIB_PREFIX}scene_implicit_function_item")
  add_library(${SCENE_IMPLICIT_FUNCTION_ITEM_LIB} SHARED
    Scene_implicit_function_item.cpp Scene_implicit_function_item.moc Color_ramp.cpp)
  target_link_libraries(${SCENE_IMPLICIT_FUNCTION_ITEM_LIB} ${SCENE_ITEM_LIB} ${QGLVIEWER_LIBRARIES} ${QT_LIBRARIES})
  set_target_properties(${SCENE_IMPLICIT_FUNCTION_ITEM_LIB} PROPERTIES DEFINE_SYMBOL scene_implicit_function_item_EXPORTS)


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

  target_link_libraries( Mesh_3 ${VIEWER_LIB} )


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
    target_link_libraries( ${plugin_name} ${SCENE_ITEM_LIB})
    # Link with CGAL
    target_link_libraries( ${plugin_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
  endmacro(polyhedron_demo_plugin)

  set(IO_IMAGE_PLUGIN_LIB "${MESH_3_LIB_PREFIX}io_image_plugin")
  polyhedron_demo_plugin(${IO_IMAGE_PLUGIN_LIB} Io_image_plugin ${imgUI_FILES})
  target_link_libraries(${IO_IMAGE_PLUGIN_LIB} ${SCENE_SEGMENTED_IMAGE_ITEM_LIB})

  if(GLEW_FOUND)
    target_link_libraries(${IO_IMAGE_PLUGIN_LIB} ${GLEW_LIBRARIES})
  endif()

  set(IO_OFF_PLUGIN_LIB "${MESH_3_LIB_PREFIX}io_off_plugin")
  polyhedron_demo_plugin(${IO_OFF_PLUGIN_LIB} Io_off_plugin)
  target_link_libraries(${IO_OFF_PLUGIN_LIB} ${SCENE_POLYHEDRON_ITEM_LIB} ${POLYGON_SOUP_LIB})

  if(GLEW_FOUND)
    set(VOLUME_PLANES_PLUGIN_LIB "${MESH_3_LIB_PREFIX}volume_planes_plugin")
    polyhedron_demo_plugin(${VOLUME_PLANES_PLUGIN_LIB} Volume_planes_plugin
	                   ${VOLUME_MOC_OUTFILES}
                           Volume_plane_intersection.cpp)                      
    target_link_libraries(${VOLUME_PLANES_PLUGIN_LIB} ${SCENE_SEGMENTED_IMAGE_ITEM_LIB} ${VIEWER_LIB})
  else()
    message(STATUS "NOTICE: The Volume_plane plugin requires GLEW and will not be compiled.")
  endif()

  set(IO_C3T3_PLUGIN_LIB "${MESH_3_LIB_PREFIX}io_c3t3_plugin")
  polyhedron_demo_plugin(${IO_C3T3_PLUGIN_LIB} Io_c3t3_plugin)
  target_link_libraries(${IO_C3T3_PLUGIN_LIB} ${SCENE_C3T3_ITEM_LIB})

  if ( EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/c3t3_rib_exporter_plugin.cpp" )
    set(C3T3_RIB_EXPORTER_PLUGIN_LIB "${MESH_3_LIB_PREFIX}c3t3_rib_exporter_plugin")
    polyhedron_demo_plugin(${C3T3_RIB_EXPORTER_PLUGIN_LIB} C3t3_rib_exporter_plugin ${ribUI_FILES})
    target_link_libraries(${C3T3_RIB_EXPORTER_PLUGIN_LIB} ${SCENE_C3T3_ITEM_LIB} ${VIEWER_LIB})
  endif()

  set(IO_IMPLICIT_FUNCTION_PLUGIN_LIB "${MESH_3_LIB_PREFIX}io_implicit_function_plugin")
  polyhedron_demo_plugin(${IO_IMPLICIT_FUNCTION_PLUGIN_LIB} Io_implicit_function_plugin ${funcUI_FILES})
  target_link_libraries(${IO_IMPLICIT_FUNCTION_PLUGIN_LIB} ${SCENE_IMPLICIT_FUNCTION_ITEM_LIB})

  set(MESH_3_PLUGIN_LIB "${MESH_3_LIB_PREFIX}mesh_3_plugin") 
  polyhedron_demo_plugin(${MESH_3_PLUGIN_LIB} Mesh_3_plugin 
                         Mesh_3_plugin_polyhedron_cgal_code.cpp
                         Mesh_3_plugin_image_cgal_code.cpp
                         Mesh_3_plugin_implicit_function_cgal_code.cpp
                         Meshing_thread.cpp
                         Scene_c3t3_item.moc
                         ${meshingUI_FILES})

  target_link_libraries(${MESH_3_PLUGIN_LIB}
                        ${SCENE_C3T3_ITEM_LIB}
                        ${SCENE_POLYHEDRON_ITEM_LIB}
                        ${SCENE_SEGMENTED_IMAGE_ITEM_LIB}
                        ${SCENE_IMPLICIT_FUNCTION_ITEM_LIB}
                        ${QGLVIEWER_LIBRARIES}
                        ${OPENGL_gl_LIBRARY}
                        ${OPENGL_glu_LIBRARY})

  set(MESH_3_OPTIMIZATION_PLUGIN_LIB "${MESH_3_LIB_PREFIX}mesh_3_optimization_plugin") 
  polyhedron_demo_plugin(${MESH_3_OPTIMIZATION_PLUGIN_LIB} Mesh_3_optimization_plugin 
                         Mesh_3_optimization_plugin_cgal_code.cpp
                         Optimizer_thread.cpp
                         Scene_c3t3_item.moc
                         ${optimUI_FILES})

  target_link_libraries(${MESH_3_OPTIMIZATION_PLUGIN_LIB}
                        ${SCENE_C3T3_ITEM_LIB}
                        ${SCENE_POLYHEDRON_ITEM_LIB}
                        ${SCENE_SEGMENTED_IMAGE_ITEM_LIB}
                        ${SCENE_IMPLICIT_FUNCTION_ITEM_LIB}
                        ${QGLVIEWER_LIBRARIES}
                        ${OPENGL_gl_LIBRARY}
                        ${OPENGL_glu_LIBRARY})


else (CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND AND Boost_FOUND)

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

  if(NOT QGLVIEWER_FOUND)
    set(MESH_3_MISSING_DEPS "QGLViewer, ${MESH_3_MISSING_DEPS}")
  endif()

  if(NOT QGLVIEWER_FOUND)
    set(MESH_3_MISSING_DEPS "QGLViewer, ${MESH_3_MISSING_DEPS}")
  endif()

  if(NOT Boost_FOUND)
    set(MESH_3_MISSING_DEPS "Boost, ${MESH_3_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${MESH_3_MISSING_DEPS}and will not be compiled.")

endif (CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND AND Boost_FOUND)



