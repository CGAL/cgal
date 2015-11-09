include( polyhedron_demo_macros )


  polyhedron_demo_plugin(pca_plugin Polyhedron_demo_pca_plugin)
  target_link_libraries(pca_plugin scene_polyhedron_item scene_basic_objects)