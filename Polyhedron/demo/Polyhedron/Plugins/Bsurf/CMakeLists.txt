include(polyhedron_demo_macros)

polyhedron_demo_plugin(locally_shortest_path_plugin Locally_shortest_path_plugin)
target_link_libraries(locally_shortest_path_plugin PUBLIC scene_edit_path_item
                                                          scene_surface_mesh_item)

