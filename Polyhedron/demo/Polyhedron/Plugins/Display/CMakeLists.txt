include(polyhedron_demo_macros)

qt6_wrap_ui(display_propertyUI_FILES Display_property.ui)
polyhedron_demo_plugin(display_property_plugin Display_property_plugin ${display_propertyUI_FILES} KEYWORDS Viewer)
target_link_libraries(display_property_plugin PUBLIC scene_surface_mesh_item
                                                     scene_points_with_normal_item
                                                     scene_color_ramp)

if(TARGET CGAL::Eigen3_support)
  qt6_wrap_ui(heat_methodUI_FILES Heat_method.ui)
  polyhedron_demo_plugin(heat_method_plugin Heat_method_plugin ${heat_methodUI_FILES} KEYWORDS Viewer)
  target_link_libraries(heat_method_plugin PUBLIC scene_surface_mesh_item
                                                  scene_selection_item
                                                  scene_color_ramp
                                                  CGAL::Eigen3_support)
endif()
