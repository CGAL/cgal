include(polyhedron_demo_macros)

polyhedron_demo_plugin(do_trees_intersect_plugin Do_trees_intersect_plugin)
target_link_libraries(do_trees_intersect_plugin PUBLIC scene_surface_mesh_item
                                                       scene_movable_sm_item)

polyhedron_demo_plugin(cut_plugin Cut_plugin)
target_link_libraries(cut_plugin PUBLIC scene_surface_mesh_item
                                        scene_basic_objects scene_color_ramp)
if(TARGET CGAL::TBB_support)
  target_link_libraries(cut_plugin PUBLIC CGAL::TBB_support)
endif()
