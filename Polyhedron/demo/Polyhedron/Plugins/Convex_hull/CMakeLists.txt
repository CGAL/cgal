include(polyhedron_demo_macros)

polyhedron_demo_plugin(convex_hull_plugin Convex_hull_plugin KEYWORDS
                       PointSetProcessing)
target_link_libraries(
  convex_hull_plugin PUBLIC scene_points_with_normal_item scene_polylines_item
                            scene_selection_item scene_surface_mesh_item)

polyhedron_demo_plugin(kernel_plugin Kernel_plugin)
target_link_libraries(kernel_plugin PUBLIC scene_surface_mesh_item)
