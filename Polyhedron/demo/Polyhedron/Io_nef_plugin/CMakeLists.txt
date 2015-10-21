include( polyhedron_demo_macros )


  polyhedron_demo_plugin(nef_io_plugin Polyhedron_demo_io_nef_plugin)
  target_link_libraries(nef_io_plugin  scene_nef_polyhedron_item)