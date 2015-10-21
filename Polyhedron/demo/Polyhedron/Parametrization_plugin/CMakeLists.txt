include( polyhedron_demo_macros )

    polyhedron_demo_plugin(parameterization_plugin Polyhedron_demo_parameterization_plugin)
    target_link_libraries(parameterization_plugin scene_polyhedron_item scene_textured_polyhedron_item )
