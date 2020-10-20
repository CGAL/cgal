include(polyhedron_demo_macros)
if(TARGET CGAL::Eigen_support)
  qt5_wrap_ui( display_propertyUI_FILES Display_property.ui )
  polyhedron_demo_plugin(display_property_plugin Display_property_plugin ${display_propertyUI_FILES} KEYWORDS Viewer)
  target_link_libraries(display_property_plugin
    PUBLIC
    scene_surface_mesh_item
    scene_points_with_normal_item
    scene_color_ramp
    CGAL::Eigen_support)
endif()
