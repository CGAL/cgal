set(list_of_whitelisted_headers_txt [=[
  CGAL/SCIP_mixed_integer_program_traits.h
  CGAL/GLPK_mixed_integer_program_traits.h
  CGAL/Classification/TensorFlow/Neural_network_classifier.h
  CGAL/Linear_cell_complex_constructors.h
  CGAL/CGAL_Ipelet_base.h
  CGAL/IO/read_las_points.h
  CGAL/IO/write_las_points.h
  CGAL/IO/read_ply_points.h
  CGAL/IO/write_ply_points.h
  CGAL/Surface_mesh_parameterization/internal/shortest_path.h
  CGAL/Three/Edge_container.h
  CGAL/Three/exceptions.h
  CGAL/Three/Point_container.h
  CGAL/Three/Polyhedron_demo_io_plugin_interface.h
  CGAL/Three/Polyhedron_demo_plugin_helper.h
  CGAL/Three/Polyhedron_demo_plugin_interface.h
  CGAL/Three/Primitive_container.h
  CGAL/Three/Scene_draw_interface.h
  CGAL/Three/Scene_group_item.h
  CGAL/Three/Scene_interface.h
  CGAL/Three/Scene_item_config.h
  CGAL/Three/Scene_item.h
  CGAL/Three/Scene_item_rendering_helper.h
  CGAL/Three/Scene_item_with_properties.h
  CGAL/Three/Scene_print_item_interface.h
  CGAL/Three/Scene_transparent_interface.h
  CGAL/Three/Scene_zoomable_item_interface.h
  CGAL/Three/TextRenderer.h
  CGAL/Three/Triangle_container.h
  CGAL/Three/Viewer_config.h
  CGAL/Three/Viewer_interface.h
  CGAL/Three/Three.h

]=])

separate_arguments(list_of_whitelisted_headers UNIX_COMMAND ${list_of_whitelisted_headers_txt})
