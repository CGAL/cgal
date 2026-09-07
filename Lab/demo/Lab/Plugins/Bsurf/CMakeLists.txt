include(CGALlab_macros)

if(TARGET CGAL::Eigen3_support)
  cgal_lab_plugin(locally_shortest_path_plugin Locally_shortest_path_plugin)
  target_link_libraries(locally_shortest_path_plugin PUBLIC locally_shortest_path_item
                                                            scene_surface_mesh_item
                                                            scene_polylines_item
                                                            CGAL::Eigen3_support)
endif()
