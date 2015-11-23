include( polyhedron_demo_macros )

 polyhedron_demo_plugin(group_plugin Polyhedron_demo_group_plugin ${groupUI_FILES})
 target_link_libraries(group_plugin scene_group_item)