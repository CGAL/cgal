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

#option(POLYHEDRON_DEMO_ENABLE_FORWARD_DECL "In the Polyhedron demo, enable " OFF)
#mark_as_advanced(POLYHEDRON_DEMO_ENABLE_FORWARD_DECL)

# Let plugins be compiled in the same directory as the executable.
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

# Include this package's headers first
include_directories( BEFORE ./ ./include ../../include ./CGAL_demo )

add_subdirectory( implicit_functions )

# Find CGAL

option( POLYHEDRON_QTSCRIPT_DEBUGGER 
  "Activate the use of Qt Script Debugger in Polyhedron_3 demo" OFF)

# Find Qt5 itself
set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )
if( POLYHEDRON_QTSCRIPT_DEBUGGER)
  set( QT_USE_QTSCRIPTTOOLS  TRUE )
endif()
find_package(CGAL COMPONENTS Qt5)
include( ${CGAL_USE_FILE} )
find_package(Qt5 QUIET COMPONENTS OpenGL Script Svg Xml ScriptTools)

# Find OpenGL
find_package(OpenGL)

# Find QGLViewer
if(Qt5_FOUND)

  add_definitions(-DQT_NO_KEYWORDS)
  include_directories( ${QT_INCLUDE_DIR}/QtScriptTools )
  find_package(QGLViewer )

endif(Qt5_FOUND)

find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
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

if(CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND)

  set(Boost_USE_MULTITHREADED      ON)
  find_package(Boost COMPONENTS thread system)

  include_directories ( ${QGLVIEWER_INCLUDE_DIR} )

  qt5_wrap_ui( MainWindowUI_files MainWindow.ui )
  qt5_wrap_ui( FileLoaderDialogUI_files FileLoaderDialog.ui )
  qt5_wrap_ui( Show_point_dialogUI_FILES Show_point_dialog.ui )
  qt5_wrap_ui( remeshingUI_FILES  Remeshing_dialog.ui)
  qt5_wrap_ui( Mean_curvature_flow_skeleton_pluginUI_FILES  Mean_curvature_flow_skeleton_plugin.ui)
  qt5_wrap_ui( meshingUI_FILES  Meshing_dialog.ui )
  qt5_wrap_ui( cameraUI_FILES   Camera_positions_list.ui )
  qt5_wrap_ui( PreferencesUI_FILES   Preferences.ui )
  qt5_wrap_ui( point_inside_polyhedronUI_FILES Point_inside_polyhedron_widget.ui)
  qt5_wrap_ui( point_set_selectionUI_FILES Point_set_selection_widget.ui)
  qt5_wrap_ui( polyhedron_slicerUI_FILES Polyhedron_slicer_widget.ui)
  qt5_wrap_ui( segmentationUI_FILES Mesh_segmentation_widget.ui)
  qt5_wrap_ui( hole_fillingUI_FILES Hole_filling_widget.ui)
  qt5_wrap_ui( fairingUI_FILES Fairing_widget.ui)
  qt5_wrap_ui( selectionUI_FILES Selection_widget.ui)
  qt5_wrap_ui( funcUI_FILES Function_dialog.ui )
  qt5_wrap_ui( shortestPathUI_FILES Shortest_path_widget.ui )

  qt5_generate_moc( "MainWindow.h" "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp" )
  qt5_generate_moc( "File_loader_dialog.h" "${CMAKE_CURRENT_BINARY_DIR}/File_loader_dialog_moc.cpp" )
  add_file_dependencies( MainWindow_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h" )
  add_file_dependencies( File_loader_dialog_moc_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/File_loader_dialog.h" )

  qt5_generate_moc( "Viewer.h" "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp" )
  add_file_dependencies( Viewer_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h" )

  include( ${CMAKE_CURRENT_SOURCE_DIR}/polyhedron_demo_macros.cmake )


  qt5_add_resources ( RESOURCE_FILES Polyhedron_3.qrc )
  qt5_add_resources(gl_splat_rc GlSplat/glsplat.qrc)
  add_library(gl_splat SHARED
      GlSplat/GlSplat.cpp GlSplat/Shader.cpp ${gl_splat_rc})
  qt5_use_modules(gl_splat OpenGL Gui Xml)
  target_link_libraries(gl_splat ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} ${QT_QTOPENGL_LIBRARY} ${GLEW_LIBRARIES})

  qt5_generate_moc( "CGAL_demo/Viewer_interface.h" "${CMAKE_CURRENT_BINARY_DIR}/Viewer_interface_moc.cpp" )


  # AUXILIARY LIBRARIES

  # put plugins (which are shared libraries) at the same location as
  # executable files
  set(CGAL_POLYHEDRON_DEMO_PLUGINS_DIR "${RUNTIME_OUTPUT_PATH}")
  set(LIBRARY_OUTPUT_PATH "${CGAL_POLYHEDRON_DEMO_PLUGINS_DIR}")

  add_library(demo_framework SHARED
    Scene.cpp
    Viewer.cpp #Viewer_moc.cpp
    Viewer_interface_moc.cpp
    Scene_item.cpp #Scene_item.moc
    Polyhedron_demo_plugin_helper.cpp)
    qt5_use_modules(demo_framework OpenGL Gui Widgets Script Xml)

  target_link_libraries(demo_framework
    ${QGLVIEWER_LIBRARIES}
    ${QT_LIBRARIES}
    ${OPENGL_gl_LIBRARY}
    ${OPENGL_glu_LIBRARY}
    )

  add_library(scene_basic_objects SHARED
    Scene_plane_item.cpp
    #Scene_plane_item.moc
)
  target_link_libraries(scene_basic_objects 
    demo_framework
    ${CGAL_LIBRARIES}
    ${QGLVIEWER_LIBRARIES}
    ${QT_LIBRARIES}
    ${OPENGL_gl_LIBRARY}
    ${OPENGL_glu_LIBRARY}
    )
  qt5_use_modules(scene_basic_objects OpenGL Gui Xml Script Widgets)
  add_library(point_dialog SHARED Show_point_dialog.cpp Show_point_dialog.ui ${Show_point_dialogUI_FILES})
  qt5_use_modules(point_dialog OpenGL Gui Xml Script Widgets)
  target_link_libraries(point_dialog ${QT_LIBRARIES})

  macro(add_item item_name)
    add_library(${item_name} SHARED ${ARGN})
    qt5_use_modules(${item_name} OpenGL Gui Xml Script Widgets)
    target_link_libraries(${item_name} demo_framework ${QT_LIBRARIES} ${CGAL_LIBRARIES} ${Boost_LIBRARIES})
  endmacro(add_item)

  add_item(scene_c2t3_item Scene_c2t3_item.cpp)# Scene_c2t3_item.moc)
  add_item(scene_polyhedron_item Scene_polyhedron_item.cpp)# Scene_polyhedron_item.moc)
  add_item(scene_polyhedron_transform_item Scene_polyhedron_transform_item.cpp )#Scene_polyhedron_transform_item.moc)
  # special
  target_link_libraries(scene_polyhedron_transform_item scene_polyhedron_item)

  add_item(scene_combinatorial_map_item Scene_combinatorial_map_item.cpp)# Scene_combinatorial_map_item.moc)
  # special
  target_link_libraries(scene_combinatorial_map_item scene_polyhedron_item)

  add_item(scene_polylines_item Scene_polylines_item.cpp)# Scene_polylines_item.moc)

  add_item(scene_polyhedron_item_decorator Scene_polyhedron_item_decorator.cpp )#Scene_polyhedron_item_decorator.moc)
  target_link_libraries(scene_polyhedron_item_decorator scene_polyhedron_item)
  
  add_item(scene_polyhedron_item_k_ring_selection Scene_polyhedron_item_k_ring_selection.cpp)# Scene_polyhedron_item_k_ring_selection.moc)
  target_link_libraries(scene_polyhedron_item_k_ring_selection scene_polyhedron_item)
 
  add_item(scene_polyhedron_selection_item Scene_polyhedron_selection_item.cpp)# Scene_polyhedron_selection_item.moc)
  target_link_libraries(scene_polyhedron_selection_item scene_polyhedron_item_decorator scene_polyhedron_item_k_ring_selection)

  add_item(scene_polyhedron_shortest_path_item Scene_polyhedron_shortest_path_item.cpp)
  
  target_link_libraries(scene_polyhedron_shortest_path_item scene_polyhedron_item_decorator scene_polyhedron_item scene_polylines_item)

  if(EIGEN3_FOUND )
    add_item(scene_textured_polyhedron_item Scene_textured_polyhedron_item.cpp texture.cpp)# Scene_textured_polyhedron_item.moc)
  endif()

  add_item(scene_implicit_function_item Scene_implicit_function_item.cpp  Color_ramp.cpp )#Scene_implicit_function_item.moc)

  add_item(scene_polygon_soup_item Scene_polygon_soup_item.cpp)# Scene_polygon_soup_item.moc)
  target_link_libraries(scene_polygon_soup_item scene_polyhedron_item)
  add_item(scene_nef_polyhedron_item Scene_nef_polyhedron_item.cpp)# Scene_nef_polyhedron_item.moc
  
  target_link_libraries(scene_nef_polyhedron_item scene_polyhedron_item)
  add_item(scene_points_with_normal_item Scene_points_with_normal_item.cpp)# Scene_points_with_normal_item.moc)
  target_link_libraries( scene_points_with_normal_item gl_splat ${GLEW_LIBRARIES} )
  target_link_libraries( demo_framework gl_splat ${GLEW_LIBRARIES} )

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
    set(CGAL_POLYHEDRON_DEMO_DEFINITIONS "-DUSE_FORWARD_DECL")
    add_definitions(${CGAL_POLYHEDRON_DEMO_DEFINITIONS})
    add_executable  ( Polyhedron_3
      MainWindow.cpp
      Polyhedron_3.cpp
      File_loader_dialog_moc.cpp
       ${FileLoaderDialogUI_files} ${MainWindowUI_files} ${PreferencesUI_FILES} ${RESOURCE_FILES} )
    add_to_cached_list( CGAL_EXECUTABLE_TARGETS Polyhedron_3 )
    qt5_use_modules(Polyhedron_3 Gui OpenGL Xml Widgets Script Svg )
    if(EIGEN3_FOUND)
#      add_executable( Polyhedron_3 Scene_tex_rendering.cpp Scene_tex_polyhedron_operations.cpp )
    endif()
#  else(POLYHEDRON_DEMO_ENABLE_FORWARD_DECL)
#    add_file_dependencies( Polyhedron_3.cpp "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp"
#                                            "${CMAKE_CURRENT_BINARY_DIR}/Scene_moc.cpp"
#                                            "${CMAKE_CURRENT_BINARY_DIR}/File_loader_dialog_moc.cpp"
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
  target_link_libraries( Polyhedron_3 gl_splat ${GLEW_LIBRARIES})

  # Link with CGAL
  target_link_libraries( Polyhedron_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  # Link with libQGLViewer, OpenGL
  target_link_libraries( Polyhedron_3 ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )


  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Polyhedron_3 )


  ###########
  # PLUGINS #
  ###########
  remove_definitions(-DQT_STATICPLUGIN)

  polyhedron_demo_plugin(remeshing_plugin Polyhedron_demo_remeshing_plugin Polyhedron_demo_remeshing_plugin_cgal_code.cpp Mesher_base.cpp ${remeshingUI_FILES})
  target_link_libraries(remeshing_plugin scene_polyhedron_item scene_polygon_soup_item scene_c2t3_item)

  polyhedron_demo_plugin(camera_positions_plugin 
    Polyhedron_demo_camera_positions_plugin
    Camera_positions_list 
    #Camera_positions_list.moc
    ${cameraUI_FILES})

  target_link_libraries( camera_positions_plugin ${QGLVIEWER_LIBRARIES} demo_framework)

  if ( Boost_VERSION GREATER 103400 )
    qt5_generate_moc( "Polyhedron_demo_mesh_3_plugin_cgal_code.cpp" "${CMAKE_CURRENT_BINARY_DIR}/Scene_c3t3_item.moc" )

    polyhedron_demo_plugin(mesh_3_plugin Polyhedron_demo_mesh_3_plugin 
      Polyhedron_demo_mesh_3_plugin_cgal_code.cpp #Scene_c3t3_item.moc
      ${meshingUI_FILES})
    target_link_libraries(mesh_3_plugin scene_polyhedron_item scene_polygon_soup_item ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})
  else( Boost_VERSION GREATER 103400 )
    message(STATUS "warning: the plugin mesh_3_plugin requires Boost>=1.34.1 and will not be compiled.")
  endif( Boost_VERSION GREATER 103400 )

  polyhedron_demo_plugin(inside_out_plugin Polyhedron_demo_inside_out_plugin)
  target_link_libraries(inside_out_plugin scene_polyhedron_item scene_polygon_soup_item)

  polyhedron_demo_plugin(off_plugin Polyhedron_demo_off_plugin)
  target_link_libraries(off_plugin scene_polyhedron_item scene_polygon_soup_item)

  polyhedron_demo_plugin(nef_io_plugin Polyhedron_demo_io_nef_plugin)
  target_link_libraries(nef_io_plugin  scene_nef_polyhedron_item)

  polyhedron_demo_plugin(polylines_io_plugin Polyhedron_demo_polylines_io_plugin)
  target_link_libraries(polylines_io_plugin scene_polylines_item)

  polyhedron_demo_plugin(io_implicit_function_plugin Io_implicit_function_plugin ${funcUI_FILES})
  target_link_libraries(io_implicit_function_plugin scene_implicit_function_item)

  polyhedron_demo_plugin(gocad_plugin Polyhedron_demo_gocad_plugin)
  target_link_libraries(gocad_plugin scene_polyhedron_item)

  polyhedron_demo_plugin(stl_plugin Polyhedron_demo_stl_plugin)
  target_link_libraries(stl_plugin scene_polyhedron_item scene_polygon_soup_item)

  polyhedron_demo_plugin(xyz_plugin Polyhedron_demo_xyz_plugin)
  target_link_libraries(xyz_plugin scene_points_with_normal_item)
  
  polyhedron_demo_plugin(selection_io_plugin Polyhedron_demo_selection_io_plugin)
  target_link_libraries(selection_io_plugin scene_polyhedron_selection_item)

  polyhedron_demo_plugin(orient_soup_plugin Polyhedron_demo_orient_soup_plugin)
  target_link_libraries(orient_soup_plugin scene_polygon_soup_item scene_polyhedron_item)

  polyhedron_demo_plugin(triangulate_facets_plugin Polyhedron_demo_triangulate_facets_plugin)
  target_link_libraries(triangulate_facets_plugin scene_polyhedron_item)

  polyhedron_demo_plugin(off_to_nef_plugin Polyhedron_demo_off_to_nef_plugin)
  target_link_libraries(off_to_nef_plugin scene_nef_polyhedron_item)

  polyhedron_demo_plugin(off_to_xyz_plugin Polyhedron_demo_off_to_xyz_plugin)
  target_link_libraries(off_to_xyz_plugin scene_points_with_normal_item)

  polyhedron_demo_plugin(convex_hull_plugin Polyhedron_demo_convex_hull_plugin)
  target_link_libraries(convex_hull_plugin scene_polyhedron_item scene_points_with_normal_item scene_polylines_item scene_polyhedron_selection_item)

  polyhedron_demo_plugin(kernel_plugin Polyhedron_demo_kernel_plugin)
  target_link_libraries(kernel_plugin scene_polyhedron_item)

  polyhedron_demo_plugin(pca_plugin Polyhedron_demo_pca_plugin)
  target_link_libraries(pca_plugin scene_polyhedron_item scene_basic_objects)

  qt5_wrap_ui( advancing_frontUI_FILES Polyhedron_demo_advancing_front_plugin.ui)
    polyhedron_demo_plugin(advancing_front_plugin Polyhedron_demo_advancing_front_plugin ${advancing_frontUI_FILES})
    target_link_libraries(advancing_front_plugin scene_polygon_soup_item scene_points_with_normal_item)

  if(EIGEN3_FOUND)
    qt5_wrap_ui( scale_spaceUI_FILES Polyhedron_demo_scale_space_reconstruction_plugin.ui)
    polyhedron_demo_plugin(scale_space_reconstruction_plugin Polyhedron_demo_scale_space_reconstruction_plugin ${scale_spaceUI_FILES})
    target_link_libraries(scale_space_reconstruction_plugin scene_polygon_soup_item scene_polyhedron_item scene_points_with_normal_item)

    polyhedron_demo_plugin(parameterization_plugin Polyhedron_demo_parameterization_plugin)
    target_link_libraries(parameterization_plugin scene_polyhedron_item scene_textured_polyhedron_item )

    qt5_wrap_ui( poissonUI_FILES Polyhedron_demo_poisson_plugin.ui)
    polyhedron_demo_plugin(poisson_plugin Polyhedron_demo_poisson_plugin Polyhedron_demo_poisson_plugin_impl.cpp ${poissonUI_FILES})
    target_link_libraries(poisson_plugin scene_polyhedron_item scene_points_with_normal_item)

    qt5_wrap_ui( normal_estimationUI_FILES Polyhedron_demo_normal_estimation_plugin.ui)
    polyhedron_demo_plugin(normal_estimation_plugin Polyhedron_demo_normal_estimation_plugin ${normal_estimationUI_FILES})
    target_link_libraries(normal_estimation_plugin scene_points_with_normal_item)

    qt5_wrap_ui( vcm_normal_estimationUI_FILES Polyhedron_demo_vcm_normal_estimation_plugin.ui)
    polyhedron_demo_plugin(vcm_normal_estimation_plugin Polyhedron_demo_vcm_normal_estimation_plugin ${vcm_normal_estimationUI_FILES})
    target_link_libraries(vcm_normal_estimation_plugin scene_points_with_normal_item)

    qt5_wrap_ui( features_detection_pluginUI_FILES Polyhedron_demo_features_detection_plugin.ui)
    polyhedron_demo_plugin(features_detection_plugin Polyhedron_demo_features_detection_plugin ${features_detection_pluginUI_FILES})
    target_link_libraries(features_detection_plugin scene_points_with_normal_item)
    
    polyhedron_demo_plugin(point_set_smoothing_plugin Polyhedron_demo_point_set_smoothing_plugin)
    target_link_libraries(point_set_smoothing_plugin scene_points_with_normal_item)

    polyhedron_demo_plugin(point_set_average_spacing_plugin Polyhedron_demo_point_set_average_spacing_plugin)
    target_link_libraries(point_set_average_spacing_plugin scene_points_with_normal_item)

    polyhedron_demo_plugin(jet_fitting_plugin Polyhedron_demo_jet_fitting_plugin)
    target_link_libraries(jet_fitting_plugin scene_polyhedron_item scene_polylines_item)
  else(EIGEN3_FOUND)
    message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Parameterization plugin will not be available.")
    message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Poisson reconstruction plugin will not be available.")
    message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Normal estimation plugins will not be available.")
    message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Smoothing plugin will not be available.")
    message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Average spacing plugin will not be available.")
    message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Jet fitting plugin will not be available.")
    message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Scale space reconstruction plugin will not be available.")
    message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Feature detection plugin will not be available.")
  endif(EIGEN3_FOUND)

  polyhedron_demo_plugin(self_intersection_plugin Polyhedron_demo_self_intersection_plugin)
  target_link_libraries(self_intersection_plugin scene_polyhedron_item scene_polyhedron_selection_item)

  polyhedron_demo_plugin(polyhedron_stitching_plugin Polyhedron_demo_polyhedron_stitching_plugin)
  target_link_libraries(polyhedron_stitching_plugin scene_polyhedron_item scene_polylines_item)

  polyhedron_demo_plugin(join_and_split_polyhedra_plugin Polyhedron_demo_join_and_split_polyhedra_plugin)
  target_link_libraries(join_and_split_polyhedra_plugin scene_polyhedron_item)

  polyhedron_demo_plugin(subdivision_methods_plugin Polyhedron_demo_subdivision_methods_plugin)
  target_link_libraries(subdivision_methods_plugin scene_polyhedron_item)

  polyhedron_demo_plugin(mesh_simplification_plugin Polyhedron_demo_mesh_simplification_plugin)
  target_link_libraries(mesh_simplification_plugin scene_polyhedron_item)

  polyhedron_demo_plugin(polyhedron_slicer_plugin Polyhedron_demo_polyhedron_slicer_plugin ${polyhedron_slicerUI_FILES})
  target_link_libraries(polyhedron_slicer_plugin scene_polyhedron_item scene_basic_objects scene_polylines_item)

  polyhedron_demo_plugin(mesh_segmentation_plugin Polyhedron_demo_mesh_segmentation_plugin ${segmentationUI_FILES})
  target_link_libraries(mesh_segmentation_plugin scene_polyhedron_item)
  
  polyhedron_demo_plugin(mean_curvature_flow_skeleton_plugin Polyhedron_demo_mean_curvature_flow_skeleton_plugin ${Mean_curvature_flow_skeleton_pluginUI_FILES})
  target_link_libraries(mean_curvature_flow_skeleton_plugin 
    scene_polyhedron_item 
    scene_points_with_normal_item
    scene_polylines_item
    demo_framework)

  polyhedron_demo_plugin(nef_plugin Polyhedron_demo_nef_plugin)
  target_link_libraries(nef_plugin scene_nef_polyhedron_item)

  polyhedron_demo_plugin(intersection_plugin Polyhedron_demo_intersection_plugin)
  target_link_libraries(intersection_plugin  scene_polyhedron_item scene_polylines_item)

  polyhedron_demo_plugin(transform_polyhedron_plugin Polyhedron_demo_transform_polyhedron_plugin)
  target_link_libraries(transform_polyhedron_plugin  scene_polyhedron_item scene_polyhedron_transform_item)

  polyhedron_demo_plugin(corefinement_plugin Polyhedron_demo_corefinement_plugin)
  target_link_libraries(corefinement_plugin  scene_polyhedron_item scene_combinatorial_map_item scene_polylines_item)

  polyhedron_demo_plugin(trivial_plugin Polyhedron_demo_trivial_plugin)

  # Edit polyhedron scene item and plugin
  if ( EIGEN3_FOUND AND "${EIGEN3_VERSION}" VERSION_GREATER "3.1.90" )
    qt5_wrap_ui( editionUI_FILES Deform_mesh.ui )
    #qt5_automoc( Scene_edit_polyhedron_item.cpp )
    add_item(scene_edit_polyhedron_item Scene_edit_polyhedron_item.cpp #Scene_edit_polyhedron_item.moc
      ${editionUI_FILES})

    target_link_libraries(scene_edit_polyhedron_item scene_polyhedron_item scene_polyhedron_item_k_ring_selection)

    polyhedron_demo_plugin(edit_polyhedron_plugin Polyhedron_demo_edit_polyhedron_plugin ${editionUI_FILES})
    target_link_libraries(edit_polyhedron_plugin scene_polyhedron_item scene_edit_polyhedron_item)
  else()
    message(STATUS "NOTICE: The polyhedron edit plugin require Eigen 3.2 (or higher) and will not be available.")
  endif()
  
  polyhedron_demo_plugin(cut_plugin Polyhedron_demo_cut_plugin)
  target_link_libraries(cut_plugin scene_polyhedron_item scene_basic_objects)
  
  polyhedron_demo_plugin(point_inside_polyhedron_plugin Polyhedron_demo_point_inside_polyhedron_plugin ${point_inside_polyhedronUI_FILES})
  target_link_libraries(point_inside_polyhedron_plugin scene_polyhedron_item scene_points_with_normal_item)
  
  qt5_wrap_ui(point_set_simplificationUI_FILES Polyhedron_demo_point_set_simplification_plugin.ui)
  polyhedron_demo_plugin(point_set_simplification_plugin Polyhedron_demo_point_set_simplification_plugin ${point_set_simplificationUI_FILES})
  target_link_libraries(point_set_simplification_plugin scene_points_with_normal_item)

  qt5_wrap_ui( ps_outliers_removal_UI_FILES Polyhedron_demo_point_set_outliers_removal_plugin.ui)

  polyhedron_demo_plugin(point_set_selection_plugin Polyhedron_demo_point_set_selection_plugin ${point_set_selectionUI_FILES})
  target_link_libraries(point_set_selection_plugin scene_points_with_normal_item scene_polylines_item)
  
  qt5_wrap_ui(point_set_shape_detectionUI_FILES Polyhedron_demo_point_set_shape_detection_plugin.ui)
  polyhedron_demo_plugin(point_set_shape_detection_plugin Polyhedron_demo_point_set_shape_detection_plugin ${point_set_shape_detectionUI_FILES})
  target_link_libraries(point_set_shape_detection_plugin scene_polyhedron_item scene_points_with_normal_item)

  qt5_wrap_ui(point_set_wlopFILES Polyhedron_demo_point_set_wlop_plugin.ui)
  polyhedron_demo_plugin(point_set_wlop_plugin Polyhedron_demo_point_set_wlop_plugin ${point_set_wlopFILES})
  target_link_libraries(point_set_wlop_plugin scene_points_with_normal_item)

  qt5_wrap_ui(point_set_bilateral_smoothingUI_FILES Polyhedron_demo_point_set_bilateral_smoothing_plugin.ui)
  polyhedron_demo_plugin(point_set_bilateral_smoothing_plugin Polyhedron_demo_point_set_bilateral_smoothing_plugin ${point_set_bilateral_smoothingUI_FILES})
  target_link_libraries(point_set_bilateral_smoothing_plugin scene_points_with_normal_item)

  qt5_wrap_ui(point_set_upsamplingUI_FILES Polyhedron_demo_point_set_upsampling_plugin.ui)
  polyhedron_demo_plugin(point_set_upsampling_plugin Polyhedron_demo_point_set_upsampling_plugin ${point_set_upsamplingUI_FILES})
  target_link_libraries(point_set_upsampling_plugin scene_points_with_normal_item)

  qt5_wrap_ui( ps_outliers_removal_UI_FILES Polyhedron_demo_point_set_outliers_removal_plugin.ui)

  polyhedron_demo_plugin(point_set_outliers_removal_plugin Polyhedron_demo_point_set_outliers_removal_plugin ${ps_outliers_removal_UI_FILES}) 
  target_link_libraries(point_set_outliers_removal_plugin scene_points_with_normal_item)
  
  if(EIGEN3_FOUND AND ${EIGEN3_VERSION} VERSION_GREATER "3.1.90")
    polyhedron_demo_plugin(hole_filling_plugin Polyhedron_demo_hole_filling_plugin ${hole_fillingUI_FILES})
    target_link_libraries(hole_filling_plugin scene_polyhedron_item scene_polylines_item)

    polyhedron_demo_plugin(hole_filling_polyline_plugin Polyhedron_demo_hole_filling_polyline_plugin)
    target_link_libraries(hole_filling_polyline_plugin scene_polyhedron_item scene_polylines_item)

    polyhedron_demo_plugin(fairing_plugin Polyhedron_demo_fairing_plugin ${fairingUI_FILES})
    target_link_libraries(fairing_plugin scene_polyhedron_selection_item)
  else()
    message(STATUS "NOTICE: The hole filling and fairing plugins require Eigen 3.2 (or higher) and will not be available.")
  endif()

  polyhedron_demo_plugin(selection_plugin Polyhedron_demo_selection_plugin ${selectionUI_FILES})
  target_link_libraries(selection_plugin scene_polyhedron_selection_item scene_points_with_normal_item scene_polylines_item)
  
  polyhedron_demo_plugin(shortest_path_plugin Polyhedron_demo_shortest_path_plugin ${shortestPathUI_FILES})
  target_link_libraries(shortest_path_plugin scene_polyhedron_item scene_polylines_item scene_polyhedron_selection_item scene_polyhedron_shortest_path_item scene_basic_objects)
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
