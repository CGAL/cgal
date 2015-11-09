include( polyhedron_demo_macros )

qt5_wrap_ui( Mean_curvature_flow_skeleton_pluginUI_FILES  Mean_curvature_flow_skeleton_plugin.ui)
  polyhedron_demo_plugin(mean_curvature_flow_skeleton_plugin Polyhedron_demo_mean_curvature_flow_skeleton_plugin ${Mean_curvature_flow_skeleton_pluginUI_FILES})
  target_link_libraries(mean_curvature_flow_skeleton_plugin 
    scene_polyhedron_item 
    scene_points_with_normal_item
    scene_polylines_item
    demo_framework)