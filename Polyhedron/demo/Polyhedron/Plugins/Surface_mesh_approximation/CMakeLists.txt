include( polyhedron_demo_macros )

qt5_wrap_ui(basicUI_FILES Surface_mesh_approximation_dockwidget.ui)
polyhedron_demo_plugin(surface_mesh_approximation_plugin
  Surface_mesh_approximation_plugin
  ${basicUI_FILES}
  Color_cheat_sheet)
target_link_libraries(surface_mesh_approximation_plugin
  PUBLIC
  scene_polyhedron_item
  scene_surface_mesh_item
  scene_polygon_soup_item
  scene_points_with_normal_item
  scene_polylines_item
  scene_polyhedron_selection_item)
