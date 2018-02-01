set(list_of_whitelisted_headers_txt [=[
  CGAL/Linear_cell_complex_constructors.h
  CGAL/CGAL_Ipelet_base.h
  CGAL/IO/read_las_points.h
  CGAL/IO/write_las_points.h
  CGAL/IO/read_ply_points.h
  CGAL/IO/write_ply_points.h
  CGAL/Surface_mesh_parameterization/internal/shortest_path.h
exceptions.h
Polyhedron_demo_plugin_interface.h
Scene_interface.h    Scene_item_with_properties.h   Scene_zoomable_item_interface.h  Viewer_interface.h
Polyhedron_demo_io_plugin_interface.h  Scene_draw_interface.h              Scene_item_config.h  Scene_print_item_interface.h   TextRenderer.h
Polyhedron_demo_plugin_helper.h        Scene_group_item.h                  Scene_item.h         Scene_transparent_interface.h  Viewer_config.h

]=])

separate_arguments(list_of_whitelisted_headers UNIX_COMMAND ${list_of_whitelisted_headers_txt})
