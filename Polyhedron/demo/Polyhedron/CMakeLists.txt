# This is the CMake script for compiling the CGAL Polyhedron demo.

project( Polyhedron )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

#option(POLYHEDRON_DEMO_ENABLE_FORWARD_DECL "In the Polyhedron demo, enable " OFF)
#mark_as_advanced(POLYHEDRON_DEMO_ENABLE_FORWARD_DECL)

# Let plugins be compiled in the same directory as the executable.
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

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

option( POLYHEDRON_QTSCRIPT_DEBUGGER 
  "Activate the use of Qt Script Debugger in Polyhedron_3 demo" OFF)

# Find Qt4 itself
set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )
if( POLYHEDRON_QTSCRIPT_DEBUGGER)
  set( QT_USE_QTSCRIPTTOOLS  TRUE )
endif()
find_package(Qt4)

# Find OpenGL
find_package(OpenGL)

# Find QGLViewer
if(QT4_FOUND)
  include(${QT_USE_FILE})
  include_directories( ${QT_INCLUDE_DIR}/QtScriptTools )
  find_package(QGLViewer )
endif(QT4_FOUND)



# Eigen is now used by default
find_package(Eigen3 3.0.91) #(requires 3.1.0-alpha1 or greater)
if (NOT EIGEN3_FOUND)
  # Find LAPACK (optional), for curvatures estimation
  find_package(LAPACK)
  if(LAPACK_FOUND)
    include( ${LAPACK_USE_FILE} )
  endif(LAPACK_FOUND)

  # Find TAUCS (optionnal), for parametrization
  find_package(TAUCS)
  if(TAUCS_FOUND)
    include( ${TAUCS_USE_FILE} )
  endif(TAUCS_FOUND)
else()
  include( ${EIGEN3_USE_FILE} )
endif()


if(CGAL_Qt4_FOUND AND QT4_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  find_package(Boost)

  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )

  # Parameterization needs Eigen3 or TAUCS
  if(NOT EIGEN3_FOUND AND NOT TAUCS_FOUND)
    message(STATUS "NOTICE: Eigen 3.1 (or greater) and TAUCS is not found. parametrization will not be available.")
  endif(NOT EIGEN3_FOUND AND NOT TAUCS_FOUND)

  # Curvature estimation needs Eigen3 or LAPACK
  if(NOT EIGEN3_FOUND AND NOT LAPACK_FOUND)
    message(STATUS "NOTICE: Eigen 3.1 (or greater) and LAPACK is not found. curvatures estimation will not be available.")
  endif(NOT EIGEN3_FOUND AND NOT LAPACK_FOUND)

  qt4_wrap_ui( MainWindowUI_files MainWindow.ui )
  qt4_wrap_ui( Show_point_dialogUI_FILES Show_point_dialog.ui )
  qt4_wrap_ui( remeshingUI_FILES  Remeshing_dialog.ui)
  qt4_wrap_ui( meshingUI_FILES  Meshing_dialog.ui Meshing_pause_widget.ui )

  include(AddFileDependencies)

  qt4_generate_moc( "MainWindow.h" "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp" )
  add_file_dependencies( MainWindow_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" )

  qt4_generate_moc( "Viewer.h" "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp" )
  add_file_dependencies( Viewer_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h" )

  qt4_generate_moc( "Polyhedron_demo_remeshing_plugin_cgal_code.cpp" "${CMAKE_CURRENT_BINARY_DIR}/Polyhedron_demo_remeshing_plugin_cgal_code.moc" )

  qt4_add_resources ( RESOURCE_FILES Polyhedron_3.qrc )

  qt4_automoc(Scene.cpp
              Scene_item.cpp
              Scene_plane_item.cpp
              Scene_polygon_soup_item.cpp 
              Scene_polyhedron_item.cpp
              Scene_edit_polyhedron_item.cpp
              Scene_textured_polyhedron_item.cpp
              Scene_c2t3_item.cpp
              Scene_nef_polyhedron_item.cpp
              Show_point_dialog.cpp
              )

  # AUXILIARY LIBRARIES

  # put plugins (which are shared libraries) at the same location as
  # executable files
  set(LIBRARY_OUTPUT_PATH ${RUNTIME_OUTPUT_PATH})

  add_library(demo_framework SHARED
    Scene.cpp
    Viewer.cpp Viewer_moc.cpp
    Scene_item.cpp Scene_item.moc
    Scene_item_with_display_list.cpp
    Polyhedron_demo_plugin_helper.cpp)

  target_link_libraries(demo_framework
    ${QGLVIEWER_LIBRARIES}
    ${OPENGL_gl_LIBRARY}
    ${OPENGL_glu_LIBRARY}
    )

  add_library(scene_basic_objects SHARED
    Scene_plane_item.cpp Scene_plane_item.moc)
  target_link_libraries(scene_basic_objects 
    demo_framework
    ${QGLVIEWER_LIBRARIES}
    ${OPENGL_gl_LIBRARY}
    ${OPENGL_glu_LIBRARY}
    )

  add_library(point_dialog SHARED Show_point_dialog.cpp Show_point_dialog.ui ${Show_point_dialogUI_FILES} )

  add_library(scene_c2t3_item SHARED
    Scene_c2t3_item.cpp Scene_c2t3_item.moc)
  target_link_libraries(scene_c2t3_item demo_framework)

  add_library(scene_polyhedron_item SHARED
    Scene_polyhedron_item.cpp Scene_polyhedron_item.moc)
  target_link_libraries(scene_polyhedron_item demo_framework)

  add_library(scene_edit_polyhedron_item SHARED
    Scene_edit_polyhedron_item.cpp Scene_edit_polyhedron_item.moc)
  target_link_libraries(scene_edit_polyhedron_item scene_polyhedron_item demo_framework)

  if(EIGEN3_FOUND OR TAUCS_FOUND)
  add_library(scene_textured_polyhedron_item SHARED
    Scene_textured_polyhedron_item.cpp texture.cpp Scene_textured_polyhedron_item.moc)
  target_link_libraries(scene_textured_polyhedron_item demo_framework)
  endif(EIGEN3_FOUND OR TAUCS_FOUND)

  add_library(scene_polygon_soup_item SHARED
    Scene_polygon_soup_item.cpp Scene_polygon_soup_item.moc)
  target_link_libraries(scene_polygon_soup_item demo_framework scene_polyhedron_item)

  add_library(scene_nef_polyhedron_item SHARED
    Scene_nef_polyhedron_item.cpp Scene_nef_polyhedron_item.moc
    Scene_nef_rendering.cpp)
  target_link_libraries(scene_nef_polyhedron_item demo_framework scene_polyhedron_item)

  foreach( lib 
      demo_framework
      scene_basic_objects
      scene_polyhedron_item
      scene_polygon_soup_item
      scene_nef_polyhedron_item )
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${lib} )
  endforeach()

  add_definitions(-DQT_STATICPLUGIN)
#  if(POLYHEDRON_DEMO_ENABLE_FORWARD_DECL)
    add_definitions(-DUSE_FORWARD_DECL)
    add_executable  ( Polyhedron_3 MainWindow.cpp
      Polyhedron_3.cpp
#      Viewer.cpp
#      Scene.cpp
#      MainWindow_curvature_estimation.cpp
      MainWindow_moc.cpp
#      Viewer_moc.cpp
      ${MainWindowUI_files} ${RESOURCE_FILES} )
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS Polyhedron_3 )
    if(EIGEN3_FOUND OR TAUCS_FOUND)
#      add_executable( Polyhedron_3 Scene_tex_rendering.cpp Scene_tex_polyhedron_operations.cpp )
    endif()
#  else(POLYHEDRON_DEMO_ENABLE_FORWARD_DECL)
#    add_file_dependencies( Polyhedron_3.cpp "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp"
#                                            "${CMAKE_CURRENT_BINARY_DIR}/Scene_moc.cpp"
#					    "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp" )
#    add_executable  ( Polyhedron_3 Polyhedron_3.cpp ${UI_FILES} ${RESOURCE_FILES} )
#  endif(POLYHEDRON_DEMO_ENABLE_FORWARD_DECL)


  # Link with Qt libraries
  if( POLYHEDRON_QTSCRIPT_DEBUGGER )
    add_definitions( -DQT_SCRIPTTOOLS_LIB )
    target_link_libraries( Polyhedron_3 QtScriptTools )
  endif()
  target_link_libraries( Polyhedron_3 ${QT_LIBRARIES} )
  target_link_libraries( Polyhedron_3 demo_framework )
  target_link_libraries( Polyhedron_3 point_dialog )

  # Link with CGAL
  target_link_libraries( Polyhedron_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  # Link with libQGLViewer, OpenGL
  target_link_libraries( Polyhedron_3 ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Polyhedron_3 )


  ###########
  # PLUGINS #
  ###########
  remove_definitions(-DQT_STATICPLUGIN)

  macro(polyhedron_demo_plugin plugin_name plugin_implementation_base_name)
    list_split(option ARGN_TAIL ${ARGN} )
    if(NOT ${option} STREQUAL "EXCLUDE_FROM_ALL")
      if(NOT ${option} STREQUAL "NO_MOC")
        set(other_sources ${ARGN})
        set(option "")
      else()
        set(other_sources ${ARGN_TAIL})
      endif()
    else()
      set(other_sources ${ARGN_TAIL})
    endif()
    if("${option}" STREQUAL "NO_MOC") 
      set(option "")
      set(moc_file_name "")
    else()
      set(moc_file_name ${plugin_implementation_base_name}.moc )
      qt4_generate_moc( ${plugin_implementation_base_name}.cpp "${CMAKE_CURRENT_BINARY_DIR}/${moc_file_name}" )
      add_file_dependencies( ${moc_file_name} "${CMAKE_CURRENT_SOURCE_DIR}/${plugin_implementation_base_name}.cpp" )
    endif()

    add_library(${plugin_name} MODULE ${option} ${moc_file_name} ${plugin_implementation_base_name}.cpp ${other_sources})
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${plugin_name} )
    # Link with Qt
    target_link_libraries( ${plugin_name} ${QT_LIBRARIES} )
    # Link with the demo_framework
    target_link_libraries( ${plugin_name} demo_framework)
    # Link with CGAL
    target_link_libraries( ${plugin_name} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
    add_dependencies( ${plugin_name} Polyhedron_3 )
  endmacro(polyhedron_demo_plugin)

  polyhedron_demo_plugin(remeshing_plugin Polyhedron_demo_remeshing_plugin Polyhedron_demo_remeshing_plugin_cgal_code.cpp Polyhedron_demo_remeshing_plugin_cgal_code.moc ${remeshingUI_FILES})
  target_link_libraries(remeshing_plugin scene_polyhedron_item scene_polygon_soup_item scene_c2t3_item)

  if ( Boost_VERSION GREATER 103400 )
    qt4_generate_moc( "Polyhedron_demo_mesh_3_plugin_cgal_code.cpp" "${CMAKE_CURRENT_BINARY_DIR}/Scene_c3t3_item.moc" )

    polyhedron_demo_plugin(mesh_3_plugin Polyhedron_demo_mesh_3_plugin 
      Polyhedron_demo_mesh_3_plugin_cgal_code.cpp Scene_c3t3_item.moc)
    target_link_libraries(mesh_3_plugin scene_polyhedron_item ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})
  else( Boost_VERSION GREATER 103400 )
    message(STATUS "warning: the plugin mesh_3_plugin requires Boost>=1.34.1 and will not be compiled.")
  endif( Boost_VERSION GREATER 103400 )

  polyhedron_demo_plugin(inside_out_plugin Polyhedron_demo_inside_out_plugin)
  target_link_libraries(inside_out_plugin scene_polyhedron_item scene_polygon_soup_item)

  polyhedron_demo_plugin(off_plugin Polyhedron_demo_off_plugin)
  target_link_libraries(off_plugin scene_polyhedron_item scene_polygon_soup_item)

  polyhedron_demo_plugin(orient_soup_plugin Polyhedron_demo_orient_soup_plugin)
  target_link_libraries(orient_soup_plugin scene_polygon_soup_item scene_polyhedron_item)

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

  if(EIGEN3_FOUND OR TAUCS_FOUND)
    polyhedron_demo_plugin(parameterization_plugin Polyhedron_demo_parameterization_plugin)
    target_link_libraries(parameterization_plugin scene_polyhedron_item scene_textured_polyhedron_item )
  endif(EIGEN3_FOUND OR TAUCS_FOUND)

  polyhedron_demo_plugin(self_intersection_plugin Polyhedron_demo_self_intersection_plugin)
  target_link_libraries(self_intersection_plugin scene_polyhedron_item)

  polyhedron_demo_plugin(subdivision_methods_plugin Polyhedron_demo_subdivision_methods_plugin)
  target_link_libraries(subdivision_methods_plugin scene_polyhedron_item)

  polyhedron_demo_plugin(mesh_simplification_plugin Polyhedron_demo_mesh_simplification_plugin)
  target_link_libraries(mesh_simplification_plugin scene_polyhedron_item)

  polyhedron_demo_plugin(nef_plugin Polyhedron_demo_nef_plugin)
  target_link_libraries(nef_plugin scene_nef_polyhedron_item)

  polyhedron_demo_plugin(trivial_plugin Polyhedron_demo_trivial_plugin)

  polyhedron_demo_plugin(edit_polyhedron_plugin Polyhedron_demo_edit_polyhedron_plugin)
  target_link_libraries(edit_polyhedron_plugin scene_polyhedron_item scene_edit_polyhedron_item)

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
