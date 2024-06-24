include(CGALlab_macros)

cgal_lab_plugin(locally_shortest_path_plugin Locally_shortest_path_plugin)
target_link_libraries(locally_shortest_path_plugin PUBLIC locally_shortest_path_item
                                                          scene_surface_mesh_item
                                                          scene_polylines_item)


