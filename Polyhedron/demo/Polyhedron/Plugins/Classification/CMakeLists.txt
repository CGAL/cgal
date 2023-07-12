include(polyhedron_demo_macros)

if(TARGET CGAL::Eigen3_support)

  find_package(Boost QUIET OPTIONAL_COMPONENTS serialization iostreams)
  include(CGAL_Boost_serialization_support)
  include(CGAL_Boost_iostreams_support)
  if(NOT TARGET CGAL::Boost_serialization_support OR NOT TARGET CGAL::Boost_iostreams_support)
    message(STATUS "NOTICE: Boost IO Streams and/or Serialization not found, reading deprecated Classification config files won't be possible.")
  endif()

  find_package(OpenCV QUIET COMPONENTS core ml) # Need core + machine learning
  set_package_properties(
    OpenCV PROPERTIES
    DESCRIPTION "A library for real-time computer vision."
    PURPOSE "Enables the random forest predicate for the classification plugin."
  )
  include(CGAL_OpenCV_support)
  if(NOT TARGET CGAL::OpenCV_support)
    message(STATUS "NOTICE: OpenCV was not found. OpenCV random forest predicate for classification won't be available.")
  endif()

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
    CGAL::Eigen3_support)

    if(BUILD_TESTING AND NOT CMAKE_VS_MSBUILD_COMMAND)
      set_tests_properties(
        "compilation of  classification_plugin"
        PROPERTIES RESOURCE_LOCK Selection_test_resources)
    endif()

  if(TARGET CGAL::Boost_serialization_support AND TARGET CGAL::Boost_iostreams_support)
    target_link_libraries(classification_plugin PUBLIC CGAL::Boost_serialization_support
                                                       CGAL::Boost_iostreams_support)
  endif()

  if(TARGET CGAL::OpenCV_support)
    target_link_libraries(classification_plugin PUBLIC CGAL::OpenCV_support)
  endif()

  if(TARGET CGAL::TBB_support)
    target_link_libraries(classification_plugin PUBLIC CGAL::TBB_support)
  endif()

  add_dependencies(classification_plugin point_set_selection_plugin selection_plugin)

else()
  message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Classification plugin will not be available.")
endif()
