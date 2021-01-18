include(polyhedron_demo_macros)

if(TARGET CGAL::Eigen_support)

  set(Classification_dependencies_met TRUE)

  find_package(Boost OPTIONAL_COMPONENTS serialization iostreams)
  include(CGAL_Boost_serialization_support)
  include(CGAL_Boost_iostreams_support)
  if(NOT TARGET CGAL::Boost_serialization_support)
    message(
      STATUS
        "NOTICE: Boost Serialization not found. Classification plugin won't be available."
    )
    set(Classification_dependencies_met FALSE)
  endif()

  if(NOT TARGET CGAL::Boost_iostreams_support)
    message(
      STATUS
        "NOTICE: Boost IOStreams not found. Classification plugin won't be available."
    )
    set(Classification_dependencies_met FALSE)
  endif()

  find_package(OpenCV QUIET COMPONENTS core ml) # Need core + machine learning
  set_package_properties(
    OpenCV PROPERTIES
    DESCRIPTION "A library for  real-time computer vision."
    PURPOSE "Enables the random forest predicate for the classification plugin."
  )
  include(CGAL_OpenCV_support)
  if(NOT TARGET CGAL::OpenCV_support)
    message(
      STATUS
        "NOTICE: OpenCV was not found. OpenCV random forest predicate for classification won't be available."
    )

  endif()

  find_package(TensorFlow QUIET)
  set_package_properties(
    TensorFlow PROPERTIES
    DESCRIPTION "A library for machine learning."
    PURPOSE
      "Enables the Neural Network predicate for the classification plugin.")

  include(CGAL_TensorFlow_support)
  if(NOT TARGET CGAL::TensorFlow_support)
    message(
      STATUS
        "NOTICE: TensorFlow not found, Neural Network predicate for classification won't be available."
    )
  endif()

  if(Classification_dependencies_met)
    qt5_wrap_ui(classificationUI_FILES Classification_widget.ui
                Classification_advanced_widget.ui)
    polyhedron_demo_plugin(
      classification_plugin
      Classification_plugin
      Point_set_item_classification.cpp
      Cluster_classification.cpp
      Surface_mesh_item_classification.cpp
      ${classificationUI_FILES}
      KEYWORDS
      Classification)
    target_link_libraries(
      classification_plugin
      PUBLIC scene_points_with_normal_item
             scene_polylines_item
             scene_polygon_soup_item
             scene_surface_mesh_item
             scene_selection_item
             scene_color_ramp
             CGAL::Eigen_support
             CGAL::Boost_serialization_support
             CGAL::Boost_iostreams_support)

    if(OpenCV_FOUND)
      target_link_libraries(classification_plugin PUBLIC CGAL::OpenCV_support)
    endif()
    if(TensorFlow_FOUND)
      target_link_libraries(classification_plugin
                            PUBLIC CGAL::TensorFlow_support)
    endif()
    if(TBB_FOUND)
      target_link_libraries(classification_plugin PUBLIC CGAL::TBB_support)
    endif()

    add_dependencies(classification_plugin point_set_selection_plugin
                     selection_plugin)
  endif()

else()
  message(
    STATUS
      "NOTICE: Eigen 3.1 (or greater) was not found. Classification plugin will not be available."
  )
endif()
