# This is the CMake script for compiling the CGAL Polyhedron demo.

project( Polyhedron )
# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

# Compatibility with CMake 3.0
if(POLICY CMP0042)
  # Do not enable the use of MACOSX_RPATH
  # http://www.cmake.org/cmake/help/v3.0/policy/CMP0042.html
  cmake_policy(SET CMP0042 OLD)
endif()

#Defines flags to emulate windows behavior for linking error generation
if(CMAKE_CXX_COMPILER_ID EQUAL Clang OR CMAKE_COMPILER_IS_GNUCC  OR CMAKE_COMPILER_IS_GNUCXX)
  if(UNIX OR APPLE)
    SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fvisibility=hidden" )
  endif()
  if(UNIX AND NOT APPLE)
    SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -z defs")
    SET( CMAKE_MODULE_LINKER_FLAGS  "${CMAKE_MODULE_LINKER_FLAGS} -z defs")
    SET( CMAKE_SHARED_LINKER_FLAGS  "${CMAKE_SHARED_LINKER_FLAGS} -z defs")
  endif()
endif()
# Let plugins be compiled in the same directory as the executable.
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

# Include this package's headers first
include_directories( BEFORE ./ ./include ../../include ./CGAL_demo )
list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_SOURCE_DIR}")

add_subdirectory( implicit_functions )

# Find CGAL

option( POLYHEDRON_QTSCRIPT_DEBUGGER 
  "Activate the use of Qt Script Debugger in Polyhedron_3 demo" OFF)

# Find Qt5 itself
find_package(CGAL COMPONENTS Qt5 ImageIO)
include( ${CGAL_USE_FILE} )

find_package(Qt5
             QUIET
             COMPONENTS            OpenGL Script Svg Xml
             OPTIONAL_COMPONENTS   ScriptTools)

# Find OpenGL
find_package(OpenGL)

# Find QGLViewer
if(Qt5_FOUND)

  add_definitions(-DQT_NO_KEYWORDS)
  add_definitions(-DSCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE)
  find_package(QGLViewer )

endif(Qt5_FOUND)

find_package(Eigen3 3.2.0) #(requires 3.2.0 or greater)
if (EIGEN3_FOUND)
  include( ${EIGEN3_USE_FILE} )
endif(EIGEN3_FOUND)

# Activate concurrency?
option(POLYHEDRON_DEMO_ACTIVATE_CONCURRENCY
  "Enable concurrency"
  ON)

if( POLYHEDRON_DEMO_ACTIVATE_CONCURRENCY )
  find_package( TBB )
  if( TBB_FOUND )
    include( ${TBB_USE_FILE} )
    list( APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES} )
  else()
    message( STATUS "NOTICE: Intel TBB was not found. Bilateral smoothing and WLOP plugins are faster if TBB is linked." )
  endif()
endif()


# Activate concurrency ? (turned OFF by default)
option(ACTIVATE_CONCURRENT_MESH_3
  "Activate parallelism in Mesh_3"
  OFF)

# And add -DCGAL_CONCURRENT_MESH_3 if that option is ON
if( ACTIVATE_CONCURRENT_MESH_3 OR ENV{ACTIVATE_CONCURRENT_MESH_3} )
  add_definitions( -DCGAL_CONCURRENT_MESH_3 )
  if(NOT TBB_FOUND)
    find_package( TBB REQUIRED )
    if( TBB_FOUND )
      include(${TBB_USE_FILE})
      list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})
    else()
      message(STATUS "NOTICE: Intel TBB was not found. Mesh_3 is faster if TBB is linked.")
    endif()
  endif()

else( ACTIVATE_CONCURRENT_MESH_3 OR ENV{ACTIVATE_CONCURRENT_MESH_3} )
  option( LINK_WITH_TBB
    "Link with TBB anyway so we can use TBB timers for profiling"
    ON)
  if( LINK_WITH_TBB )
    find_package( TBB )
  endif( LINK_WITH_TBB )
endif()




if(CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  set(Boost_USE_MULTITHREADED      ON)
  find_package(Boost COMPONENTS thread system)

  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )
  qt5_wrap_ui( MainWindowUI_files MainWindow.ui)
  qt5_wrap_ui( statisticsUI_FILES Statistics_on_item_dialog.ui)
  qt5_wrap_ui( FileLoaderDialogUI_files FileLoaderDialog.ui )
  qt5_wrap_ui( Show_point_dialogUI_FILES Show_point_dialog.ui )
  qt5_wrap_ui( PreferencesUI_FILES   Preferences.ui )
  qt5_wrap_ui( Show_point_dialogUI_FILES Show_point_dialog.ui )
  qt5_generate_moc( "File_loader_dialog.h" "${CMAKE_CURRENT_BINARY_DIR}/File_loader_dialog_moc.cpp" )
  add_file_dependencies( File_loader_dialog_moc_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/File_loader_dialog.h" )

  qt5_generate_moc( "Viewer.h" "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp" )
  add_file_dependencies( Viewer_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h" )

  include( ${CMAKE_CURRENT_SOURCE_DIR}/polyhedron_demo_macros.cmake )


  qt5_add_resources ( RESOURCE_FILES Polyhedron_3.qrc )
  qt5_add_resources(gl_splat_rc GlSplat/glsplat.qrc)
  add_library(gl_splat SHARED
      GlSplat/GlSplat.cpp GlSplat/Shader.cpp ${gl_splat_rc})
  qt5_use_modules(gl_splat OpenGL Gui Xml)

  find_path(CGAL_THREE_VIEWER_INTERFACE_H_PATH
    NAME CGAL/Three/Viewer_interface.h
    PATHS ${CGAL_INCLUDE_DIRS}
    NO_DEFAULT_PATH
    DOC "Path to CGAL/Three/Viewer_interface.h")
  hide_variable(CGAL_THREE_VIEWER_INTERFACE_H_PATH)
  if(CGAL_THREE_VIEWER_INTERFACE_H_PATH)
    qt5_generate_moc( "${CGAL_THREE_VIEWER_INTERFACE_H_PATH}/CGAL/Three/Viewer_interface.h"
      "${CMAKE_CURRENT_BINARY_DIR}/Viewer_interface_moc.cpp" )
  else()
    message(FATAL_ERROR "Cannot find <CGAL/Three/Viewer_interface.h>")
  endif()

  find_path(CGAL_THREE_SCENE_ITEM_H_PATH
    NAME CGAL/Three/Scene_item.h
    PATHS ${CGAL_INCLUDE_DIRS}
    NO_DEFAULT_PATH
    DOC "Path to CGAL/Three/Scene_item.h")
  hide_variable(CGAL_THREE_SCENE_ITEM_H_PATH)

  find_path(CGAL_THREE_SCENE_GROUP_ITEM_H_PATH
    NAME CGAL/Three/Scene_group_item.h
    PATHS ${CGAL_INCLUDE_DIRS}
    NO_DEFAULT_PATH
    DOC "Path to CGAL/Three/Scene_group_item.h")
  hide_variable(CGAL_THREE_SCENE_GROUP_ITEM_H_PATH)

  if(CGAL_THREE_SCENE_ITEM_H_PATH)
    qt5_generate_moc( "${CGAL_THREE_SCENE_ITEM_H_PATH}/CGAL/Three/Scene_item.h"
      "${CMAKE_CURRENT_BINARY_DIR}/Scene_item_moc.cpp" )
  else()
    message(FATAL_ERROR "Cannot find <CGAL/Three/Viewer_interface.h>")
  endif()

  if(CGAL_THREE_SCENE_GROUP_ITEM_H_PATH)
    qt5_generate_moc( "${CGAL_THREE_SCENE_GROUP_ITEM_H_PATH}/CGAL/Three/Scene_group_item.h"
      "${CMAKE_CURRENT_BINARY_DIR}/Scene_group_item_moc.cpp" )
  else()
    message(FATAL_ERROR "Cannot find <CGAL/Three/Scene_group_item.h>")
  endif()


  # AUXILIARY LIBRARIES

  # put s (which are shared libraries) at the same location as
  # executable files
  set(CGAL_POLYHEDRON_DEMO_PLUGINS_DIR "${RUNTIME_OUTPUT_PATH}")
  set(LIBRARY_OUTPUT_PATH "${CGAL_POLYHEDRON_DEMO_PLUGINS_DIR}")

  add_library(demo_framework SHARED
    Scene.cpp
    Viewer.cpp
    Viewer_interface_moc.cpp
    Scene_item_moc.cpp
    Scene_item.cpp
    Scene_group_item.cpp
    Scene_group_item_moc.cpp
    Polyhedron_demo_plugin_helper.cpp)
    qt5_use_modules(demo_framework OpenGL Gui Widgets Script Xml)

  target_link_libraries(demo_framework
    ${QGLVIEWER_LIBRARIES}
    ${OPENGL_gl_LIBRARY}
    ${OPENGL_glu_LIBRARY}
    )

  add_library(scene_basic_objects SHARED
    Scene_plane_item.cpp
)
  target_link_libraries(scene_basic_objects 
    demo_framework
    ${CGAL_LIBRARIES}
    ${QGLVIEWER_LIBRARIES}
    ${OPENGL_gl_LIBRARY}
    ${OPENGL_glu_LIBRARY}
    )
  qt5_use_modules(scene_basic_objects OpenGL Gui Xml Script Widgets)
  add_library(point_dialog SHARED Show_point_dialog.cpp Show_point_dialog.ui ${Show_point_dialogUI_FILES})
  qt5_use_modules(point_dialog OpenGL Gui Xml Script Widgets)

  macro(add_item item_name)
    add_library(${item_name} SHARED ${ARGN})
    qt5_use_modules(${item_name} OpenGL Gui Xml Script Widgets)
    target_link_libraries(${item_name} demo_framework ${CGAL_LIBRARIES} ${Boost_LIBRARIES})
  endmacro(add_item)

  add_item(scene_c2t3_item Scene_c2t3_item.cpp)
  add_item(scene_c3t3_item Scene_c3t3_item.cpp)
  target_link_libraries(scene_c3t3_item scene_polyhedron_item scene_polygon_soup_item ${TBB_LIBRARIES})
  add_item(scene_polyhedron_item Scene_polyhedron_item.cpp)
  add_item(scene_polyhedron_transform_item Plugins/PCA/Scene_polyhedron_transform_item.cpp )

  add_item(scene_segmented_image_item Scene_segmented_image_item.cpp)
 
  # special
  target_link_libraries(scene_polyhedron_transform_item scene_polyhedron_item)

  add_item(scene_combinatorial_map_item Plugins/Operations_on_polyhedra/Scene_combinatorial_map_item.cpp)
  # special
  target_link_libraries(scene_combinatorial_map_item scene_polyhedron_item)

  add_item(scene_polylines_item Scene_polylines_item.cpp)

  add_item(scene_polyhedron_item_decorator Scene_polyhedron_item_decorator.cpp )
  target_link_libraries(scene_polyhedron_item_decorator scene_polyhedron_item)
  
  add_item(scene_polyhedron_item_k_ring_selection Scene_polyhedron_item_k_ring_selection.cpp)
  target_link_libraries(scene_polyhedron_item_k_ring_selection scene_polyhedron_item)
 
  add_item(scene_polyhedron_selection_item Scene_polyhedron_selection_item.cpp)
  target_link_libraries(scene_polyhedron_selection_item scene_polyhedron_item_decorator scene_polyhedron_item_k_ring_selection)

  add_item(scene_polyhedron_shortest_path_item Plugins/Surface_mesh/Scene_polyhedron_shortest_path_item.cpp)
  
  target_link_libraries(scene_polyhedron_shortest_path_item scene_polyhedron_item_decorator scene_polyhedron_item scene_polylines_item)

  if(EIGEN3_FOUND )
    qt5_wrap_ui( editionUI_FILES Plugins/Surface_mesh_deformation/Deform_mesh.ui )
    add_item(scene_textured_polyhedron_item Scene_textured_polyhedron_item.cpp texture.cpp)
    add_item(scene_edit_polyhedron_item Plugins/Surface_mesh_deformation/Scene_edit_polyhedron_item.cpp
           ${editionUI_FILES})
    target_link_libraries(scene_edit_polyhedron_item scene_polyhedron_item scene_polyhedron_item_k_ring_selection)
  endif()

  add_item(scene_implicit_function_item Scene_implicit_function_item.cpp  Color_ramp.cpp )

  add_item(scene_polygon_soup_item Scene_polygon_soup_item.cpp)
  target_link_libraries(scene_polygon_soup_item scene_polyhedron_item)
  add_item(scene_nef_polyhedron_item Scene_nef_polyhedron_item.cpp)
  
  target_link_libraries(scene_nef_polyhedron_item scene_polyhedron_item)
  add_item(scene_points_with_normal_item Scene_points_with_normal_item.cpp)
  target_link_libraries( scene_points_with_normal_item gl_splat)
  target_link_libraries( demo_framework gl_splat)
  
  


  foreach( lib 
      demo_framework
      scene_basic_objects
      scene_polyhedron_item
      scene_polygon_soup_item
      scene_nef_polyhedron_item)
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${lib} )
  endforeach()

  add_definitions( -DUSE_FORWARD_DECL)
  add_library(polyhedron_demo SHARED
    MainWindow.cpp
    Polyhedron_demo.cpp
    File_loader_dialog_moc.cpp
    ${FileLoaderDialogUI_files} ${MainWindowUI_files} ${PreferencesUI_FILES} ${RESOURCE_FILES} ${statisticsUI_FILES})
  target_link_libraries(polyhedron_demo demo_framework point_dialog)
  qt5_use_modules(polyhedron_demo Gui OpenGL Xml Widgets Script Svg )
  add_executable  ( Polyhedron_3 Polyhedron_3.cpp )
  target_link_libraries( Polyhedron_3 polyhedron_demo )
  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Polyhedron_3 )

  if( POLYHEDRON_QTSCRIPT_DEBUGGER )
    if(TARGET Qt5::ScriptTools)
      qt5_use_modules(Polyhedron_3 ScriptTools)
    else()
      message(STATUS "POLYHEDRON_QTSCRIPT_DEBUGGER is set to TRUE but the Qt5 ScriptTools library was not found.")
    endif()
  endif()
  target_link_libraries( Polyhedron_3 demo_framework )

  # Link with CGAL
  target_link_libraries( Polyhedron_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  # Link with libQGLViewer, OpenGL
  target_link_libraries( Polyhedron_3 ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )


  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Polyhedron_3 )


  ###########
  # PLUGINS #
  ###########


  file(GLOB DEMO_PLUGINS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/Plugins/" "${CMAKE_CURRENT_SOURCE_DIR}/Plugins/*")
  FOREACH(SUB_DIR ${DEMO_PLUGINS})
      add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Plugins/${SUB_DIR}")
  ENDFOREACH()

#
# Exporting
#
  if(TARGET CGAL_Qt5)
    export(TARGETS CGAL CGAL_Qt5 FILE polyhedron_demo_targets.cmake NAMESPACE Polyhedron_)
  else()
    export(TARGETS FILE polyhedron_demo_targets.cmake)
  endif()
  export(
    TARGETS
    demo_framework
    scene_polyhedron_item
    scene_points_with_normal_item
    scene_implicit_function_item
    scene_polylines_item
    NAMESPACE Polyhedron_
    APPEND FILE polyhedron_demo_targets.cmake)

  export(
    TARGETS
    gl_splat
    NAMESPACE Polyhedron_
    APPEND FILE polyhedron_demo_targets.cmake)

  configure_file(CGAL_polyhedron_demoConfig.cmake.in  CGAL_polyhedron_demoConfig.cmake)
#TO DO script the activation of all the plugins.

else (CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  set(POLYHEDRON_MISSING_DEPS "")

  if(NOT CGAL_Qt5_FOUND)
    set(POLYHEDRON_MISSING_DEPS "the CGAL Qt5 library, ${POLYHEDRON_MISSING_DEPS}")
  endif()

  if(NOT Qt5_FOUND)
    set(POLYHEDRON_MISSING_DEPS "Qt5, ${POLYHEDRON_MISSING_DEPS}")
  endif()

  if(NOT OPENGL_FOUND)
    set(POLYHEDRON_MISSING_DEPS "OpenGL, ${POLYHEDRON_MISSING_DEPS}")
  endif()

  if(NOT QGLVIEWER_FOUND)
    set(POLYHEDRON_MISSING_DEPS "QGLViewer, ${POLYHEDRON_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${POLYHEDRON_MISSING_DEPS}and will not be compiled.")


endif (CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)
