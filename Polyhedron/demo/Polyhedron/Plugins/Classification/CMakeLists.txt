include( polyhedron_demo_macros )

if(EIGEN3_FOUND)

  set(Classification_dependencies_met TRUE)

  find_package( Boost OPTIONAL_COMPONENTS serialization iostreams )
  if (NOT Boost_SERIALIZATION_FOUND)
    message(STATUS "NOTICE: Boost Serialization not found. Classification plugin won't be available.")
    set(Classification_dependencies_met FALSE)
  endif()

  if (NOT Boost_IOSTREAMS_FOUND)
    message(STATUS "NOTICE: Boost IOStreams not found. Classification plugin won't be available.")
    set(Classification_dependencies_met FALSE)
  endif()

  find_package(OpenCV QUIET COMPONENTS core ml) # Need core + machine learning
  if (NOT OpenCV_FOUND)
    message(STATUS "NOTICE: OpenCV was not found. OpenCV random forest predicate for classification won't be available.")
  endif()

  find_package(TensorFlow QUIET)
  if (NOT TensorFlow_FOUND)
    message(STATUS "NOTICE: TensorFlow not found, Neural Network predicate for classification won't be available.")
  endif()

  if (Classification_dependencies_met)
    qt5_wrap_ui( classificationUI_FILES Classification_widget.ui Classification_advanced_widget.ui )
    polyhedron_demo_plugin(classification_plugin Classification_plugin Point_set_item_classification Cluster_classification Surface_mesh_item_classification ${classificationUI_FILES} KEYWORDS Classification)
    target_link_libraries(classification_plugin PUBLIC scene_points_with_normal_item
      scene_polylines_item scene_polygon_soup_item scene_surface_mesh_item scene_selection_item scene_color_ramp)

    CGAL_target_use_Eigen(classification_plugin)
    CGAL_target_use_Boost_IOStreams(classification_plugin)
    CGAL_target_use_Boost_Serialization(classification_plugin)
    if(OpenCV_FOUND)
      CGAL_target_use_OpenCV(classification_plugin)
    endif()
    if(TensorFlow_FOUND)
      CGAL_target_use_TensorFlow(classification_plugin)
    endif()
    if(TBB_FOUND)
      CGAL_target_use_TBB(classification_plugin)
    endif()

    add_dependencies(classification_plugin point_set_selection_plugin selection_plugin)
  endif()

else(EIGEN3_FOUND)
  message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Classification plugin will not be available.")
endif()


