# This is the CMake script for compiling the CGAL Mesh_3 demo.

project( Mesh_3_demo )

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


# Creates a new CMake option, turned ON by default
option(ACTIVATE_MSVC_PRECOMPILED_HEADERS
  "Activate precompiled headers in MSVC"
  OFF)

# Macro to add precompiled headers for MSVC
# This function does two things:
# 1. Enable precompiled headers on each file which is listed in "SourcesVar".
# 2. Add the content of "PrecompiledSource" (e.g. "StdAfx.cpp") to "SourcesVar".
MACRO(ADD_MSVC_PRECOMPILED_HEADER PrecompiledHeader PrecompiledSource SourcesVar)
  IF(MSVC AND ACTIVATE_MSVC_PRECOMPILED_HEADERS)
    GET_FILENAME_COMPONENT(PrecompiledBasename ${PrecompiledHeader} NAME_WE)
    SET(Sources ${${SourcesVar}})

    SET_SOURCE_FILES_PROPERTIES(${PrecompiledSource}
                                PROPERTIES COMPILE_FLAGS "/Yc\"${PrecompiledHeader}\"")
    SET_SOURCE_FILES_PROPERTIES(${Sources}
                                PROPERTIES COMPILE_FLAGS "/Yu\"${PrecompiledHeaders}\" /FI\"${PrecompiledHeader}\"")
    # Add precompiled header to SourcesVar
    LIST(APPEND ${SourcesVar} ${PrecompiledSource})
  ENDIF(MSVC AND ACTIVATE_MSVC_PRECOMPILED_HEADERS)
ENDMACRO(ADD_MSVC_PRECOMPILED_HEADER)
# The compiler might need more memory because of precompiled headers
if(MSVC AND ACTIVATE_MSVC_PRECOMPILED_HEADERS AND NOT(MSVC_VERSION LESS 1310))
  set(CGAL_C_FLAGS "${CGAL_C_FLAGS} /Zm1000")
  set(CGAL_CXX_FLAGS "${CGAL_CXX_FLAGS} /Zm1000")
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
set(CMAKE_MODULE_PATH 
  ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

# Activate concurrency ? (turned OFF by default)
option(ACTIVATE_CONCURRENT_MESH_3
  "Activate parallelism in Mesh_3"
  OFF)
  
# Find CGAL and CGAL Qt5
find_package(CGAL COMPONENTS Qt5 ImageIO)
include( ${CGAL_USE_FILE} )
# And add -DCGAL_CONCURRENT_MESH_3 if that option is ON
if( ACTIVATE_CONCURRENT_MESH_3 )
  add_definitions( -DCGAL_CONCURRENT_MESH_3 )
  find_package( TBB REQUIRED )
else( ACTIVATE_CONCURRENT_MESH_3 )
  option( LINK_WITH_TBB
    "Link with TBB anyway so we can use TBB timers for profiling"
    ON)
  if( LINK_WITH_TBB )
    find_package( TBB )
  endif( LINK_WITH_TBB )
endif()

if( TBB_FOUND )
  include(${TBB_USE_FILE})
endif()

# Find Qt5 itself
find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL)

# Find OpenGL
find_package(OpenGL)

# Find everybodies darling
find_package(Boost COMPONENTS thread system)

# Find QGLViewer

if(Qt5_FOUND)

 add_definitions(-DQT_NO_KEYWORDS)
 add_definitions(-DSCENE_SEGMENTED_IMAGE_GL_BUFFERS_AVAILABLE)

 find_package(QGLViewer )
endif(Qt5_FOUND)


if(CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND AND Boost_FOUND)

  # Add directory containing implicit function plugin source files
  add_subdirectory("implicit_functions")

  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )

  qt5_wrap_ui( UI_FILES ui_files/MainWindow.ui )
  qt5_wrap_ui( meshingUI_FILES  ui_files/Meshing_dialog.ui)
  qt5_wrap_ui( optimUI_FILES  ui_files/Smoother_dialog.ui ui_files/LocalOptim_dialog.ui)
  qt5_wrap_ui( ribUI_FILES  ui_files/Rib_dialog.ui)
  qt5_wrap_ui( funcUI_FILES  ui_files/Function_dialog.ui)
  qt5_wrap_ui( imgUI_FILES ui_files/Image_res_dialog.ui)

  include(AddFileDependencies)

  # Shared demo files directory
  set(DEMO_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/CGAL_demo")
  set(DEMO_H_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include/CGAL_demo")

  qt5_generate_moc( "MainWindow.h" "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp" )
  add_file_dependencies( MainWindow_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" )

  qt5_generate_moc( "${DEMO_H_DIR}/Viewer.h" Viewer_moc.cpp )
  add_file_dependencies( Viewer_moc.cpp "${DEMO_H_DIR}/Viewer.h" )

  qt5_generate_moc( "${DEMO_H_DIR}/Scene.h" Scene_moc.cpp )
  add_file_dependencies( Scene_moc.cpp "${DEMO_H_DIR}/Scene.h" )

  qt5_generate_moc( "${DEMO_H_DIR}/Scene_item.h" Scene_item_moc.cpp )
  add_file_dependencies( Scene_item_moc.cpp "${DEMO_H_DIR}/Scene_item.h" )

  qt5_wrap_cpp( VOLUME_MOC_OUTFILES ${CMAKE_CURRENT_SOURCE_DIR}/Volume_plane_thread.h )
  qt5_wrap_cpp( VOLUME_MOC_OUTFILES ${CMAKE_CURRENT_SOURCE_DIR}/Volume_plane_interface.h )

  qt5_add_resources ( RESOURCE_FILES Mesh_3.qrc )

  # put plugins (which are shared libraries) at the same location as
  # executable files
  set(LIBRARY_OUTPUT_PATH ${RUNTIME_OUTPUT_PATH})

  # prefix for Mesh_3 libs to avoid conflicts with Polyhedron demo
  set(MESH_3_LIB_PREFIX "mesh_3_demo_")

  # AUXILIARY LIBRARIES
  set(VIEWER_LIB "${MESH_3_LIB_PREFIX}viewer")
  add_library(${VIEWER_LIB} SHARED
    ${DEMO_SRC_DIR}/Viewer.cpp Viewer_moc.cpp)
  qt5_use_modules(${VIEWER_LIB} Xml Script OpenGL Svg)
  target_link_libraries(${VIEWER_LIB} ${QGLVIEWER_LIBRARIES} ${QT_LIBRARIES} ${OPENGL_LIBRARIES})
  set_target_properties(${VIEWER_LIB} PROPERTIES DEFINE_SYMBOL viewer_EXPORTS)

  set(SCENE_ITEM_LIB "${MESH_3_LIB_PREFIX}scene_item")
  add_library(${SCENE_ITEM_LIB} SHARED
    ${DEMO_SRC_DIR}/Scene_item.cpp
    Scene_item_moc.cpp
    ${DEMO_SRC_DIR}/Plugin_helper.cpp)
  qt5_use_modules(${SCENE_ITEM_LIB} Xml Script OpenGL Svg)
  set_target_properties(${SCENE_ITEM_LIB} PROPERTIES DEFINE_SYMBOL scene_item_EXPORTS)
  target_link_libraries(${SCENE_ITEM_LIB} ${QT_LIBRARIES} ${OPENGL_LIBRARIES})

  set(SCENE_SEGMENTED_IMAGE_ITEM_LIB "${MESH_3_LIB_PREFIX}scene_segmented_image_item")
  add_library(${SCENE_SEGMENTED_IMAGE_ITEM_LIB} SHARED
    Scene_segmented_image_item.cpp)# Scene_segmented_image_item.moc)
  qt5_use_modules(${SCENE_SEGMENTED_IMAGE_ITEM_LIB} Xml Script OpenGL Svg)
  target_link_libraries(${SCENE_SEGMENTED_IMAGE_ITEM_LIB} ${SCENE_ITEM_LIB} ${QGLVIEWER_LIBRARIES} ${CGAL_LIBRARIES} ${TBB_LIBRARIES})
  set_target_properties(${SCENE_SEGMENTED_IMAGE_ITEM_LIB} PROPERTIES DEFINE_SYMBOL scene_segmented_image_item_EXPORTS)

  set(SCENE_POLYHEDRON_ITEM_LIB "${MESH_3_LIB_PREFIX}scene_polyhedron_item")
  add_library(${SCENE_POLYHEDRON_ITEM_LIB} SHARED
    Scene_polyhedron_item.cpp )#Scene_polyhedron_item.moc)
  qt5_use_modules(${SCENE_POLYHEDRON_ITEM_LIB} Xml Script OpenGL Svg)
  target_link_libraries(${SCENE_POLYHEDRON_ITEM_LIB} ${SCENE_ITEM_LIB}  ${QGLVIEWER_LIBRARIES} ${CGAL_LIBRARIES} ${TBB_LIBRARIES})
  set_target_properties(${SCENE_POLYHEDRON_ITEM_LIB} PROPERTIES DEFINE_SYMBOL scene_polyhedron_item_EXPORTS)

  set(POLYGON_SOUP_LIB "${MESH_3_LIB_PREFIX}polygon_soup")
  add_library(${POLYGON_SOUP_LIB} SHARED
    Scene_polygon_soup.cpp )#Scene_polygon_soup.moc)
  qt5_use_modules(${POLYGON_SOUP_LIB} Xml Script OpenGL Svg)
  target_link_libraries(${POLYGON_SOUP_LIB} ${SCENE_ITEM_LIB}  ${QGLVIEWER_LIBRARIES} ${CGAL_LIBRARIES} ${TBB_LIBRARIES})
  set_target_properties(${POLYGON_SOUP_LIB} PROPERTIES DEFINE_SYMBOL polygon_soup_EXPORTS)

  set(SCENE_C3T3_ITEM_LIB "${MESH_3_LIB_PREFIX}scene_c3t3_item")
  set (SCENE_C3T3_ITEM_LIB_SOURCE_FILES Scene_c3t3_item.cpp)
  ADD_MSVC_PRECOMPILED_HEADER("StdAfx.h" "StdAfx.cpp" SCENE_C3T3_ITEM_LIB_SOURCE_FILES)
 # LIST(APPEND SCENE_C3T3_ITEM_LIB_SOURCE_FILES Scene_c3t3_item.moc)
  add_library(${SCENE_C3T3_ITEM_LIB} SHARED
    ${SCENE_C3T3_ITEM_LIB_SOURCE_FILES})
  qt5_use_modules(${SCENE_C3T3_ITEM_LIB} Xml Script OpenGL Svg)
  target_link_libraries(${SCENE_C3T3_ITEM_LIB} ${SCENE_ITEM_LIB} ${QGLVIEWER_LIBRARIES} ${QT_LIBRARIES} ${CGAL_LIBRARIES} ${Boost_LIBRARIES} ${TBB_LIBRARIES})
  set_target_properties(${SCENE_C3T3_ITEM_LIB} PROPERTIES DEFINE_SYMBOL scene_c3t3_item_EXPORTS)

  set(SCENE_IMPLICIT_FUNCTION_ITEM_LIB "${MESH_3_LIB_PREFIX}scene_implicit_function_item")
  add_library(${SCENE_IMPLICIT_FUNCTION_ITEM_LIB} SHARED
    Scene_implicit_function_item.cpp Color_ramp.cpp)#Scene_implicit_function_item.moc)
  qt5_use_modules(${SCENE_IMPLICIT_FUNCTION_ITEM_LIB} Xml Script OpenGL Svg)
  target_link_libraries(${SCENE_IMPLICIT_FUNCTION_ITEM_LIB} ${SCENE_ITEM_LIB} ${QGLVIEWER_LIBRARIES} ${QT_LIBRARIES} ${TBB_LIBRARIES})
  set_target_properties(${SCENE_IMPLICIT_FUNCTION_ITEM_LIB} PROPERTIES DEFINE_SYMBOL scene_implicit_function_item_EXPORTS)

  
  add_definitions(-DUSE_FORWARD_DECL)
  add_definitions(-DQT_STATICPLUGIN)
    set (MESH_3_SOURCE_FILES 
      MainWindow.cpp
      Mesh_3.cpp
      ${DEMO_SRC_DIR}/Scene.cpp
    #  MainWindow_moc.cpp
      Scene_moc.cpp)
    #ADD_MSVC_PRECOMPILED_HEADER("StdAfx.h" "StdAfx.cpp" MESH_3_SOURCE_FILES)
    LIST(APPEND MESH_3_SOURCE_FILES ${UI_FILES} ${RESOURCE_FILES})
    add_executable( Mesh_3 ${MESH_3_SOURCE_FILES} )
    qt5_use_modules(Mesh_3 Xml Script OpenGL Svg)
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS Mesh_3 )

  # Link with Qt libraries
  target_link_libraries( Mesh_3 ${QT_LIBRARIES} )

  # Link with CGAL
  target_link_libraries( Mesh_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${TBB_LIBRARIES})

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
    qt5_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/${plugin_implementation_base_name}.cpp" ${plugin_implementation_base_name}.moc )
    add_file_dependencies( ${plugin_implementation_base_name}.moc "${CMAKE_CURRENT_SOURCE_DIR}/${plugin_implementation_base_name}.cpp" )

    add_library(${plugin_name} MODULE ${option} ${plugin_implementation_base_name}.moc ${plugin_implementation_base_name}.cpp ${other_sources})
    qt5_use_modules(${plugin_name} Xml Script OpenGL Svg)
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${plugin_name} )
    # Link with Qt
    target_link_libraries( ${plugin_name} ${QT_LIBRARIES} )
    # Link with scene_item
    target_link_libraries( ${plugin_name} ${SCENE_ITEM_LIB})
    # Link with CGAL
    target_link_libraries( ${plugin_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${TBB_LIBRARIES})
  endmacro(polyhedron_demo_plugin)

  set(IO_IMAGE_PLUGIN_LIB "${MESH_3_LIB_PREFIX}io_image_plugin")
  polyhedron_demo_plugin(${IO_IMAGE_PLUGIN_LIB} Io_image_plugin ${imgUI_FILES})
  target_link_libraries(${IO_IMAGE_PLUGIN_LIB} ${SCENE_SEGMENTED_IMAGE_ITEM_LIB})

  set(IO_OFF_PLUGIN_LIB "${MESH_3_LIB_PREFIX}io_off_plugin")
  polyhedron_demo_plugin(${IO_OFF_PLUGIN_LIB} Io_off_plugin)
  target_link_libraries(${IO_OFF_PLUGIN_LIB} ${SCENE_POLYHEDRON_ITEM_LIB} ${POLYGON_SOUP_LIB})

  if(Qt5_FOUND)
    set(VOLUME_PLANES_PLUGIN_LIB "${MESH_3_LIB_PREFIX}volume_planes_plugin")
    polyhedron_demo_plugin(${VOLUME_PLANES_PLUGIN_LIB} Volume_planes_plugin
                           ${VOLUME_MOC_OUTFILES}
                           Volume_plane_intersection.cpp)
    target_link_libraries(${VOLUME_PLANES_PLUGIN_LIB} ${SCENE_SEGMENTED_IMAGE_ITEM_LIB} ${VIEWER_LIB})
  endif()

  set(IO_C3T3_PLUGIN_LIB "${MESH_3_LIB_PREFIX}io_c3t3_plugin")
  polyhedron_demo_plugin(${IO_C3T3_PLUGIN_LIB} Io_c3t3_plugin)
  target_link_libraries(${IO_C3T3_PLUGIN_LIB} ${SCENE_C3T3_ITEM_LIB})

  if ( EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/C3t3_rib_exporter_plugin.cpp" )
    set(C3T3_RIB_EXPORTER_PLUGIN_LIB "${MESH_3_LIB_PREFIX}c3t3_rib_exporter_plugin")
    polyhedron_demo_plugin(${C3T3_RIB_EXPORTER_PLUGIN_LIB} C3t3_rib_exporter_plugin ${ribUI_FILES})
    target_link_libraries(${C3T3_RIB_EXPORTER_PLUGIN_LIB} ${SCENE_C3T3_ITEM_LIB} ${VIEWER_LIB})
  endif()

  set(IO_IMPLICIT_FUNCTION_PLUGIN_LIB "${MESH_3_LIB_PREFIX}io_implicit_function_plugin")
  polyhedron_demo_plugin(${IO_IMPLICIT_FUNCTION_PLUGIN_LIB} Io_implicit_function_plugin ${funcUI_FILES})
  target_link_libraries(${IO_IMPLICIT_FUNCTION_PLUGIN_LIB} ${SCENE_IMPLICIT_FUNCTION_ITEM_LIB})

  set(MESH_3_PLUGIN_LIB "${MESH_3_LIB_PREFIX}mesh_3_plugin") 
  set(MESH_3_PLUGIN_SOURCE_FILES
    Mesh_3_plugin.cpp
    Mesh_3_plugin_polyhedron_cgal_code.cpp
    Mesh_3_plugin_image_cgal_code.cpp
    Mesh_3_plugin_implicit_function_cgal_code.cpp
    Meshing_thread.cpp)
  ADD_MSVC_PRECOMPILED_HEADER("StdAfx.h" "StdAfx.cpp" MESH_3_PLUGIN_SOURCE_FILES)
  LIST(REMOVE_AT MESH_3_PLUGIN_SOURCE_FILES 0) # Remove Mesh_3_plugin.cpp since it's added by polyhedron_demo_plugin
  LIST(APPEND MESH_3_PLUGIN_SOURCE_FILES  ${meshingUI_FILES})
  polyhedron_demo_plugin(${MESH_3_PLUGIN_LIB} Mesh_3_plugin ${MESH_3_PLUGIN_SOURCE_FILES})

#  set(MESH_3_PLUGIN_LIB "${MESH_3_LIB_PREFIX}mesh_3_plugin") 
#  polyhedron_demo_plugin(${MESH_3_PLUGIN_LIB} Mesh_3_plugin 
#                         Mesh_3_plugin_polyhedron_cgal_code.cpp
#                         Mesh_3_plugin_image_cgal_code.cpp
#                         Mesh_3_plugin_implicit_function_cgal_code.cpp
#                         Meshing_thread.cpp
#                         Scene_c3t3_item.moc
#                         ${meshingUI_FILES})

  target_link_libraries(${MESH_3_PLUGIN_LIB}
                        ${SCENE_C3T3_ITEM_LIB}
                        ${SCENE_POLYHEDRON_ITEM_LIB}
                        ${SCENE_SEGMENTED_IMAGE_ITEM_LIB}
                        ${SCENE_IMPLICIT_FUNCTION_ITEM_LIB}
                        ${QGLVIEWER_LIBRARIES}
                        ${OPENGL_gl_LIBRARY}
                        ${OPENGL_glu_LIBRARY}
                        ${TBB_LIBRARIES})

  set(MESH_3_OPTIMIZATION_PLUGIN_LIB "${MESH_3_LIB_PREFIX}mesh_3_optimization_plugin")
  set(MESH_3_OPTIMIZATION_PLUGIN_SOURCE_FILES
    Mesh_3_optimization_plugin.cpp
    Mesh_3_optimization_plugin_cgal_code.cpp
    Optimizer_thread.cpp)
  ADD_MSVC_PRECOMPILED_HEADER("StdAfx.h" "StdAfx.cpp" MESH_3_OPTIMIZATION_PLUGIN_SOURCE_FILES)
  LIST(REMOVE_AT MESH_3_OPTIMIZATION_PLUGIN_SOURCE_FILES 0) # Remove Mesh_3_optimization_plugin.cpp since it's added by polyhedron_demo_plugin
  LIST(APPEND MESH_3_OPTIMIZATION_PLUGIN_SOURCE_FILES ${optimUI_FILES})
  polyhedron_demo_plugin(${MESH_3_OPTIMIZATION_PLUGIN_LIB} Mesh_3_optimization_plugin ${MESH_3_OPTIMIZATION_PLUGIN_SOURCE_FILES})

  target_link_libraries(${MESH_3_OPTIMIZATION_PLUGIN_LIB}
                        ${SCENE_C3T3_ITEM_LIB}
                        ${SCENE_POLYHEDRON_ITEM_LIB}
                        ${SCENE_SEGMENTED_IMAGE_ITEM_LIB}
                        ${SCENE_IMPLICIT_FUNCTION_ITEM_LIB}
                        ${QGLVIEWER_LIBRARIES}
                        ${OPENGL_gl_LIBRARY}
                        ${OPENGL_glu_LIBRARY}
                        ${TBB_LIBRARIES})


else (CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND AND Boost_FOUND)

  set(MESH_3_MISSING_DEPS "")

  if(NOT CGAL_Qt5_FOUND)
    set(MESH_3_MISSING_DEPS "the CGAL Qt5 library, ${MESH_3_MISSING_DEPS}")
  endif()

  if(NOT Qt5_FOUND)
    set(MESH_3_MISSING_DEPS "Qt5, ${MESH_3_MISSING_DEPS}")
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

endif (CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND AND Boost_FOUND)



