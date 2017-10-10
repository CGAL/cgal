include( polyhedron_demo_macros )

if(EIGEN3_FOUND)

  qt5_wrap_ui( classificationUI_FILES Classification_widget.ui)
  polyhedron_demo_plugin(classification_plugin Classification_plugin Point_set_item_classification ${classificationUI_FILES})

  find_package(Boost COMPONENTS serialization iostreams)

  if (Boost_SERIALIZATION_FOUND)
  if (Boost_IOSTREAMS_FOUND)
  if (LASLIB_FOUND)

    find_package(OpenCV QUIET)
    if (OpenCV_FOUND)

      message(STATUS "Found OpenCV ${OpenCV_VERSION}")
      include_directories(${OpenCV_INCLUDE_DIRS})
      target_link_libraries(classification_plugin scene_points_with_normal_item scene_polylines_item scene_polygon_soup_item scene_color_ramp ${OpenCV_LIBS}
        ${Boost_SERIALIZATION_LIBRARY} ${Boost_IOSTREAMS_LIBRARY}
        ${LASLIB_LIBRARIES})

      target_compile_definitions(classification_plugin PUBLIC "-DCGAL_LINKED_WITH_OPENCV")

    else(OpenCV_FOUND)
      target_link_libraries(classification_plugin scene_points_with_normal_item
          scene_polylines_item scene_polygon_soup_item scene_color_ramp
          ${Boost_SERIALIZATION_LIBRARY} ${Boost_IOSTREAMS_LIBRARY}
          ${LASLIB_LIBRARIES})

      message(STATUS "NOTICE: OpenCV was not found. Random forest predicate for classification will not be available.")
    endif(OpenCV_FOUND)
 
  else(LASLIB_FOUND)
    message(STATUS "NOTICE: LAS was not found. Classification plugin will not be available.")
  endif(LASLIB_FOUND)
 
  else(Boost_IOSTREAMS_FOUND)
    message(STATUS "NOTICE: Boost_iostreams was not found. Classification plugin will not be  available.")
  endif(Boost_IOSTREAMS_FOUND)
 
  else(Boost_SERIALIZATION_FOUND)
    message(STATUS "NOTICE: Boost_serialization was not found. Classification plugin will not be  available.")
  endif(Boost_SERIALIZATION_FOUND)

else(EIGEN3_FOUND)
  message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Classification plugin will not be available.")
endif()

  
