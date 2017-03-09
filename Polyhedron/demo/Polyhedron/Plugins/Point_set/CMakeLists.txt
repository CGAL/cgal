include( polyhedron_demo_macros )
if(EIGEN3_FOUND)
  qt5_wrap_ui( surface_reconstructionUI_FILES Surface_reconstruction_plugin.ui)
  polyhedron_demo_plugin(surface_reconstruction_plugin Surface_reconstruction_plugin Surface_reconstruction_plugin_impl ${surface_reconstructionUI_FILES})
  target_link_libraries(surface_reconstruction_plugin scene_polygon_soup_item scene_polyhedron_item scene_points_with_normal_item)

  qt5_wrap_ui( point_set_normal_estimationUI_FILES Point_set_normal_estimation_plugin.ui)
  polyhedron_demo_plugin(point_set_normal_estimation_plugin Point_set_normal_estimation_plugin ${point_set_normal_estimationUI_FILES})
  target_link_libraries(point_set_normal_estimation_plugin scene_points_with_normal_item)

  qt5_wrap_ui( features_detection_pluginUI_FILES Features_detection_plugin.ui)
  polyhedron_demo_plugin(features_detection_plugin Features_detection_plugin ${features_detection_pluginUI_FILES})
  target_link_libraries(features_detection_plugin scene_points_with_normal_item)

  polyhedron_demo_plugin(point_set_smoothing_plugin Point_set_smoothing_plugin)
  target_link_libraries(point_set_smoothing_plugin scene_points_with_normal_item)

  polyhedron_demo_plugin(point_set_average_spacing_plugin Point_set_average_spacing_plugin)
  target_link_libraries(point_set_average_spacing_plugin scene_points_with_normal_item)

else(EIGEN3_FOUND)
  message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Surface reconstruction plugin will not be available.")
  message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Normal estimation plugins will not be available.")
  message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Smoothing plugin will not be available.")
  message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Average spacing plugin will not be available.")
  message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Feature detection plugin will not be available.")

endif()

  qt5_wrap_ui(point_set_bilateral_smoothingUI_FILES Point_set_bilateral_smoothing_plugin.ui)
  polyhedron_demo_plugin(point_set_bilateral_smoothing_plugin Point_set_bilateral_smoothing_plugin ${point_set_bilateral_smoothingUI_FILES})
  target_link_libraries(point_set_bilateral_smoothing_plugin scene_points_with_normal_item)

  qt5_wrap_ui( ps_outliers_removal_UI_FILES Point_set_outliers_removal_plugin.ui)
  polyhedron_demo_plugin(point_set_outliers_removal_plugin Point_set_outliers_removal_plugin ${ps_outliers_removal_UI_FILES})
  target_link_libraries(point_set_outliers_removal_plugin scene_points_with_normal_item)

  qt5_wrap_ui( point_set_selectionUI_FILES Point_set_selection_widget.ui)
  polyhedron_demo_plugin(point_set_selection_plugin Point_set_selection_plugin ${point_set_selectionUI_FILES})
  target_link_libraries(point_set_selection_plugin scene_points_with_normal_item scene_polylines_item scene_edit_box_item)

  qt5_wrap_ui(point_set_shape_detectionUI_FILES Point_set_shape_detection_plugin.ui)
  polyhedron_demo_plugin(point_set_shape_detection_plugin Point_set_shape_detection_plugin ${point_set_shape_detectionUI_FILES})
  target_link_libraries(point_set_shape_detection_plugin scene_surface_mesh_item scene_polyhedron_item scene_points_with_normal_item scene_polygon_soup_item)

  qt5_wrap_ui(point_set_simplificationUI_FILES Point_set_simplification_plugin.ui)
  polyhedron_demo_plugin(point_set_simplification_plugin Point_set_simplification_plugin ${point_set_simplificationUI_FILES})
  target_link_libraries(point_set_simplification_plugin scene_points_with_normal_item)


  qt5_wrap_ui(point_set_upsamplingUI_FILES Point_set_upsampling_plugin.ui)
  polyhedron_demo_plugin(point_set_upsampling_plugin Point_set_upsampling_plugin ${point_set_upsamplingUI_FILES})
  target_link_libraries(point_set_upsampling_plugin scene_points_with_normal_item)

  qt5_wrap_ui(point_set_wlopFILES Point_set_wlop_plugin.ui)
  polyhedron_demo_plugin(point_set_wlop_plugin Point_set_wlop_plugin ${point_set_wlopFILES})
  target_link_libraries(point_set_wlop_plugin scene_points_with_normal_item)

  polyhedron_demo_plugin(merge_point_sets_plugin Merge_point_sets_plugin)
  target_link_libraries(merge_point_sets_plugin scene_points_with_normal_item)

  polyhedron_demo_plugin(point_set_interference_plugin Point_set_interference_plugin)
  target_link_libraries(point_set_interference_plugin scene_points_with_normal_item)

  qt5_wrap_ui( alpha_shapeUI_FILES Alpha_shape_widget.ui )
  polyhedron_demo_plugin(alpha_shape_plugin Alpha_shape_plugin ${alpha_shapeUI_FILES})
  target_link_libraries(alpha_shape_plugin scene_points_with_normal_item scene_c3t3_item)
