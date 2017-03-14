include( polyhedron_demo_macros )

if(EIGEN3_FOUND)
  qt5_wrap_ui( classificationUI_FILES Classification_widget.ui)
  polyhedron_demo_plugin(classification_plugin Classification_plugin ${classificationUI_FILES})
  target_link_libraries(classification_plugin scene_point_set_classification_item scene_points_with_normal_item scene_polylines_item scene_polygon_soup_item)
else(EIGEN3_FOUND)
  message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Classification plugin will not be available.")
endif()

