include(CGALlab_macros)

cgal_lab_plugin(cdt_3_plugin CDT_3_plugin KEYWORDS CDT_3)
target_link_libraries(cdt_3_plugin PRIVATE scene_surface_mesh_item
                                           scene_polygon_soup_item
                                           scene_c3t3_item)
CGAL_Lab_target_use_c3t3_type(cdt_3_plugin)
