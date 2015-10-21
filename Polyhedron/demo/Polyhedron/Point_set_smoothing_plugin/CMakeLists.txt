include( polyhedron_demo_macros )


    polyhedron_demo_plugin(point_set_smoothing_plugin Polyhedron_demo_point_set_smoothing_plugin)
    target_link_libraries(point_set_smoothing_plugin scene_points_with_normal_item)