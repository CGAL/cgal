include( polyhedron_demo_macros )

if(EIGEN3_FOUND)

  qt5_wrap_ui( classificationUI_FILES Classification_widget.ui)
  polyhedron_demo_plugin(classification_plugin Classification_plugin Point_set_item_classification ${classificationUI_FILES})
  find_package(OpenCV)
  if(OpenCV_FOUND)
    target_link_libraries(classification_plugin scene_points_with_normal_item scene_polylines_item scene_polygon_soup_item scene_color_ramp ${OpenCV_LIBS})
    target_compile_definitions(classification_plugin PUBLIC "-DCGAL_LINKED_WITH_OPENCV")
  else()
    target_link_libraries(classification_plugin scene_points_with_normal_item scene_polylines_item scene_polygon_soup_item scene_color_ramp)
    message(STATUS "NOTICE: OpenCV was not found. Random forest predicate for classification will not be available.")
  endif()
  
else(EIGEN3_FOUND)
  message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Classification plugin will not be available.")
endif()

