include(CGALlab_macros)

#if the plugin has a UI file
qt6_wrap_ui(alpha_wrap_3UI_FILES alpha_wrap_3_dialog.ui)
cgal_lab_plugin(alpha_wrap_3_plugin Alpha_wrap_3_plugin ${alpha_wrap_3UI_FILES})

#if the plugin uses external libraries like scene_items
target_link_libraries(alpha_wrap_3_plugin PUBLIC scene_surface_mesh_item scene_polygon_soup_item scene_points_with_normal_item scene_selection_item scene_polylines_item)
