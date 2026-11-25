include(CGALlab_macros)

if(TARGET CGAL::Eigen3_support)
  cgal_lab_plugin(variational_medial_axis_skeleton_plugin Variational_medial_axis_skeleton_plugin)
  target_link_libraries(variational_medial_axis_skeleton_plugin PUBLIC variational_medial_axis_skeleton_item
                                                                       scene_surface_mesh_item
                                                                       CGAL::Eigen3_support)
endif()
