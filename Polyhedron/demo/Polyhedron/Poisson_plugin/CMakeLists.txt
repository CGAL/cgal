include( polyhedron_demo_macros )
 qt5_wrap_ui( poissonUI_FILES Polyhedron_demo_poisson_plugin.ui)
    polyhedron_demo_plugin(poisson_plugin Polyhedron_demo_poisson_plugin Polyhedron_demo_poisson_plugin_impl.cpp ${poissonUI_FILES})
    target_link_libraries(poisson_plugin scene_polyhedron_item scene_points_with_normal_item)
