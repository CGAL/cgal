include( polyhedron_demo_macros )

polyhedron_demo_plugin(cut_plugin Cut_plugin )
target_link_libraries(cut_plugin PUBLIC scene_polyhedron_item scene_surface_mesh_item scene_basic_objects scene_color_ramp)
if(TBB_FOUND)
  CGAL_target_use_TBB(cut_plugin)
endif()
