include( polyhedron_demo_macros )

polyhedron_demo_plugin(nef_plugin Nef_plugin CATEGORY-PLUGINS_mesh_processing)
  target_link_libraries(nef_plugin scene_nef_polyhedron_item)
