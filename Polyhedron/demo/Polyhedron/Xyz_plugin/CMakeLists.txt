include( polyhedron_demo_macros )

qt5_wrap_ui( xyzUI_FILES Add_point_set_dialog.ui)
polyhedron_demo_plugin(xyz_plugin Polyhedron_demo_xyz_plugin ${xyzUI_FILES})
target_link_libraries(xyz_plugin scene_points_with_normal_item)
