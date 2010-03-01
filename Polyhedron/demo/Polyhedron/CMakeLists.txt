# This is the CMake script for compiling the CGAL Polyhedron demo.

project( Polyhedron )

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

#option(POLYHEDRON_DEMO_ENABLE_FORWARD_DECL "In the Polyhedron demo, enable " OFF)
#mark_as_advanced(POLYHEDRON_DEMO_ENABLE_FORWARD_DECL)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# Use packages improved since CGAL 3.4
foreach(DEP_PKG AABB_tree STL_Extension GraphicsView Surface_mesher Filtered_kernel Profiling_tools Mesh_3 Triangulation_3 )
  foreach(CGAL_SVN_TRUNK ../../../ ../../trunk ..)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${CGAL_SVN_TRUNK}/${DEP_PKG}/include")
      include_directories (BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/${CGAL_SVN_TRUNK}/${DEP_PKG}/include")
      if(COMMAND break)
        break()
      endif()
    endif()
  endforeach()
endforeach()

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

# Find TAUCS (optionnal), for parametrization
find_package(TAUCS)

if(CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  find_package(Boost)

  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )

  # Link with BLAS, LAPACK and TAUCS (optional)
  if(TAUCS_FOUND)
    include( ${TAUCS_USE_FILE} )
    add_definitions(-DCGAL_TAUCS_ENABLED)
  else(TAUCS_FOUND)
    message(STATUS "NOTICE: TAUCS is not found. parametrization will not be available.")
  endif(TAUCS_FOUND)

  # Link with BLAS and LAPACK only (optional)
  if(LAPACK_FOUND)
    include( ${LAPACK_USE_FILE} )
    add_definitions(-DCGAL_LAPACK_ENABLED)
  else(LAPACK_FOUND)
    message(STATUS "NOTICE: LAPACK is not found. curvatures estimation will not be available.")
  endif(LAPACK_FOUND)

  qt4_wrap_ui( UI_FILES MainWindow.ui )

  qt4_wrap_ui( remeshingUI_FILES  Remeshing_dialog.ui)

  include(AddFileDependencies)

  qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" MainWindow_moc.cpp )
  add_file_dependencies( MainWindow_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" )

  qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h" Viewer_moc.cpp )
  add_file_dependencies( Viewer_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h" )

  qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/Scene.h" Scene_moc.cpp )
  add_file_dependencies( Scene_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Scene.h" )

  
  qt4_generate_moc( "${CMAKE_CURRENT_SOURCE_DIR}/Polyhedron_demo_remeshing_plugin_cgal_code.cpp" Polyhedron_demo_remeshing_plugin_cgal_code.moc )

  qt4_add_resources ( RESOURCE_FILES Polyhedron_3.qrc )

  qt4_automoc(Scene_item.cpp
              Scene_plane_item.cpp
              Scene_polygon_soup.cpp 
              Scene_polyhedron_item.cpp
              Scene_textured_polyhedron_item.cpp
              Scene_c2t3_item.cpp
              Scene_nef_polyhedron_item.cpp)

  # AUXILIARY LIBRARIES

  # put plugins (which are shared libraries) at the same location as
  # executable files

  set(LIBRARY_OUTPUT_PATH ${RUNTIME_OUTPUT_PATH})

  add_library(scene_item SHARED
    Scene_item.cpp Scene_item.moc
    Scene_item_with_display_list.cpp
    Polyhedron_demo_plugin_helper.cpp)

  add_library(scene_basic_objects SHARED
    Scene_plane_item.cpp Scene_plane_item.moc)
  target_link_libraries(scene_basic_objects scene_item ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})

  add_library(scene_c2t3_item SHARED
    Scene_c2t3_item.cpp Scene_c2t3_item.moc)
  target_link_libraries(scene_c2t3_item scene_item)

  add_library(scene_polyhedron_item SHARED
    Scene_polyhedron_item.cpp Scene_polyhedron_item.moc)
  target_link_libraries(scene_polyhedron_item scene_item)

  if(TAUCS_FOUND)
  add_library(scene_textured_polyhedron_item SHARED
    Scene_textured_polyhedron_item.cpp Scene_textured_polyhedron_item.moc)
  target_link_libraries(scene_textured_polyhedron_item scene_item)
  endif(TAUCS_FOUND)

  add_library(polygon_soup SHARED
    Scene_polygon_soup.cpp Scene_polygon_soup.moc)
  target_link_libraries(polygon_soup scene_item)

  add_library(scene_nef_polyhedron_item SHARED
    Scene_nef_polyhedron_item.cpp Scene_nef_polyhedron_item.moc
    Scene_nef_rendering.cpp)
  target_link_libraries(scene_nef_polyhedron_item scene_polyhedron_item)

  foreach( lib scene_item scene_basic_objects scene_polyhedron_item polygon_soup scene_nef_polyhedron_item )
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${lib} )
  endforeach()

  add_definitions(-DQT_STATICPLUGIN)
#  if(POLYHEDRON_DEMO_ENABLE_FORWARD_DECL)
    add_definitions(-DUSE_FORWARD_DECL)
    add_executable  ( Polyhedron_3 MainWindow.cpp
      Polyhedron_3.cpp
      Viewer.cpp
      Scene.cpp
#      MainWindow_curvature_estimation.cpp
      MainWindow_moc.cpp
      Scene_moc.cpp
      Viewer_moc.cpp
      texture.cpp
      ${UI_FILES} ${RESOURCE_FILES} )
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS Polyhedron_3 )
    if(TAUCS_FOUND)
#      add_executable( Polyhedron_3 Scene_tex_rendering.cpp Scene_tex_polyhedron_operations.cpp )
    endif()
#  else(POLYHEDRON_DEMO_ENABLE_FORWARD_DECL)
#    add_file_dependencies( Polyhedron_3.cpp "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp"
#                                            "${CMAKE_CURRENT_BINARY_DIR}/Scene_moc.cpp"
#					    "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp" )
#    add_executable  ( Polyhedron_3 Polyhedron_3.cpp ${UI_FILES} ${RESOURCE_FILES} )
#  endif(POLYHEDRON_DEMO_ENABLE_FORWARD_DECL)


  # Link with Qt libraries
  target_link_libraries( Polyhedron_3 ${QT_LIBRARIES} )

  # Link with CGAL
  target_link_libraries( Polyhedron_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  # Link with libQGLViewer, OpenGL
  target_link_libraries( Polyhedron_3 ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

  # Link with the scene_item library.
#  target_link_libraries( Polyhedron_3 scene_item )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Polyhedron_3 )


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
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${plugin_name} )
    # Link with Qt
    target_link_libraries( ${plugin_name} ${QT_LIBRARIES} )
    # Link with scene_item
    target_link_libraries( ${plugin_name} scene_item)
    # Link with CGAL
    target_link_libraries( ${plugin_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
  endmacro(polyhedron_demo_plugin)

  polyhedron_demo_plugin(remeshing_plugin Polyhedron_demo_remeshing_plugin Polyhedron_demo_remeshing_plugin_cgal_code.cpp Polyhedron_demo_remeshing_plugin_cgal_code.moc ${remeshingUI_FILES})
  target_link_libraries(remeshing_plugin scene_polyhedron_item polygon_soup scene_c2t3_item)

  if ( Boost_VERSION GREATER 103400 )
    qt4_generate_moc("${CMAKE_CURRENT_SOURCE_DIR}/Polyhedron_demo_mesh_3_plugin_cgal_code.cpp" Scene_c3t3_item.moc )

    polyhedron_demo_plugin(mesh_3_plugin Polyhedron_demo_mesh_3_plugin 
      Polyhedron_demo_mesh_3_plugin_cgal_code.cpp Scene_c3t3_item.moc)
    target_link_libraries(mesh_3_plugin scene_polyhedron_item ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})
  else( Boost_VERSION GREATER 103400 )
    message(STATUS "warning: the plugin mesh_3_plugin requires Boost>=1.34.1 and will not be compiled.")
  endif( Boost_VERSION GREATER 103400 )

  polyhedron_demo_plugin(inside_out_plugin Polyhedron_demo_inside_out_plugin)
  target_link_libraries(inside_out_plugin scene_polyhedron_item polygon_soup)

  polyhedron_demo_plugin(off_plugin Polyhedron_demo_off_plugin)
  target_link_libraries(off_plugin scene_polyhedron_item polygon_soup)

  polyhedron_demo_plugin(orient_soup_plugin Polyhedron_demo_orient_soup_plugin)
  target_link_libraries(orient_soup_plugin polygon_soup)

  polyhedron_demo_plugin(triangulate_facets_plugin Polyhedron_demo_triangulate_facets_plugin)
  target_link_libraries(triangulate_facets_plugin scene_polyhedron_item)

  polyhedron_demo_plugin(off_to_nef_plugin Polyhedron_demo_off_to_nef_plugin EXCLUDE_FROM_ALL)
  target_link_libraries(off_to_nef_plugin scene_nef_polyhedron_item)

  polyhedron_demo_plugin(convex_hull_plugin Polyhedron_demo_convex_hull_plugin)
  target_link_libraries(convex_hull_plugin scene_polyhedron_item)

  polyhedron_demo_plugin(kernel_plugin Polyhedron_demo_kernel_plugin)
  target_link_libraries(kernel_plugin scene_polyhedron_item)

  polyhedron_demo_plugin(pca_plugin Polyhedron_demo_pca_plugin)
  target_link_libraries(pca_plugin scene_polyhedron_item scene_basic_objects)

  if(TAUCS_FOUND)
    polyhedron_demo_plugin(parameterization_plugin Polyhedron_demo_parameterization_plugin)
    target_link_libraries(parameterization_plugin scene_polyhedron_item scene_textured_polyhedron_item )
  endif(TAUCS_FOUND)

  polyhedron_demo_plugin(self_intersection_plugin Polyhedron_demo_self_intersection_plugin)
  target_link_libraries(self_intersection_plugin scene_polyhedron_item)

  polyhedron_demo_plugin(subdivision_methods_plugin Polyhedron_demo_subdivision_methods_plugin)
  target_link_libraries(subdivision_methods_plugin scene_polyhedron_item)

  polyhedron_demo_plugin(mesh_simplification_plugin Polyhedron_demo_mesh_simplification_plugin)
  target_link_libraries(mesh_simplification_plugin scene_polyhedron_item)

  polyhedron_demo_plugin(nef_plugin Polyhedron_demo_nef_plugin)
  target_link_libraries(nef_plugin scene_nef_polyhedron_item)

  polyhedron_demo_plugin(trivial_plugin Polyhedron_demo_trivial_plugin)

  polyhedron_demo_plugin(cut_plugin Polyhedron_demo_cut_plugin)
  target_link_libraries(cut_plugin scene_polyhedron_item scene_basic_objects)

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
