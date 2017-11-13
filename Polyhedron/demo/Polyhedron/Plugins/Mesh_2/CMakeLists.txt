include( polyhedron_demo_macros )
#if the plugin has a UI file
qt5_wrap_ui( mesh_2UI_FILES mesh_2_dialog.ui)
polyhedron_demo_plugin(mesh_2_plugin Mesh_2_plugin ${mesh_2UI_FILES})
#if the plugin uses external libraries like scene_items
target_link_libraries(mesh_2_plugin scene_polyhedron_item
                                    scene_surface_mesh_item
                                    scene_polylines_item
                                    scene_points_with_normal_item)
