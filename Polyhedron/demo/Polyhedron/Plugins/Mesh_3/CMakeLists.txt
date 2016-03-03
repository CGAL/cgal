include( polyhedron_demo_macros )

if ( Boost_VERSION GREATER 103400 )
  remove_definitions(-DQT_STATICPLUGIN)

    qt5_wrap_cpp( VOLUME_MOC_OUTFILES ${CMAKE_CURRENT_SOURCE_DIR}/Volume_plane_thread.h )
    qt5_wrap_cpp( VOLUME_MOC_OUTFILES ${CMAKE_CURRENT_SOURCE_DIR}/Volume_plane_interface.h )
    qt5_wrap_ui( meshingUI_FILES  Meshing_dialog.ui Smoother_dialog.ui Local_optimizers_dialog.ui )
    qt5_generate_moc( "Polyhedron_demo_mesh_3_plugin_cgal_code.cpp" "${CMAKE_CURRENT_BINARY_DIR}/Scene_c3t3_item.moc" )
    polyhedron_demo_plugin(mesh_3_plugin Mesh_3_plugin 
      Mesh_3_plugin_cgal_code.cpp Meshing_thread.cpp 
      ${meshingUI_FILES})
    target_link_libraries(mesh_3_plugin scene_polyhedron_item scene_polygon_soup_item scene_polylines_item scene_implicit_function_item scene_segmented_image_item
	scene_c3t3_item ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})
    polyhedron_demo_plugin(mesh_3_volume_planes_plugin Volume_planes_plugin ${VOLUME_MOC_OUTFILES} Volume_plane_intersection.cpp)
    target_link_libraries(mesh_3_volume_planes_plugin scene_segmented_image_item)
    polyhedron_demo_plugin(mesh_3_optimization_plugin Optimization_plugin
      Optimization_plugin_cgal_code.cpp Optimizer_thread.cpp
      ${meshingUI_FILES})
    target_link_libraries(mesh_3_optimization_plugin scene_c3t3_item scene_polyhedron_item scene_segmented_image_item scene_implicit_function_item )

  qt5_wrap_ui( imgUI_FILES Image_res_dialog.ui raw_image.ui)
  polyhedron_demo_plugin(io_image_plugin Io_image_plugin Raw_image_dialog.cpp ${imgUI_FILES})
  target_link_libraries(io_image_plugin scene_segmented_image_item)
  polyhedron_demo_plugin(c3t3_io_plugin C3t3_io_plugin)
  target_link_libraries(c3t3_io_plugin  scene_c3t3_item)
  else( Boost_VERSION GREATER 103400 )
    message(STATUS "warning: the plugin mesh_3_plugin requires Boost>=1.34.1 and will not be compiled.")
  endif( Boost_VERSION GREATER 103400 )

  polyhedron_demo_plugin(detect_sharp_edges_plugin Detect_sharp_edges_plugin)
  target_link_libraries(detect_sharp_edges_plugin scene_polyhedron_item)
