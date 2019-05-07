include( polyhedron_demo_macros )

if(EIGEN3_FOUND)

  set(classification_linked_libraries)
  set(classification_compile_definitions)
  
  find_package(Boost OPTIONAL_COMPONENTS serialization iostreams)

  if( WIN32 )
# to avoid a warning with old cmake
    set(_Boost_BZIP2_HEADERS             "boost/iostreams/filter/bzip2.hpp")
    set(_Boost_ZLIB_HEADERS              "boost/iostreams/filter/zlib.hpp")
    find_package( Boost OPTIONAL_COMPONENTS serialization iostreams bzip2 zlib)
    if( Boost_ZLIB_FOUND AND Boost_BZIP2_FOUND)
      if(TARGET Boost::bzip2 AND TARGET Boost::zlib)
        set(classification_linked_libraries ${classification_linked_libraries}
          Boost::bzip2 Boost::zlib)
      else()
        set(classification_linked_libraries ${classification_linked_libraries}
          ${Boost_ZLIB_LIBRARY} ${Boost_BZIP2_LIBRARY})
      endif()
    endif()
  endif()

  if (Boost_SERIALIZATION_FOUND AND Boost_IOSTREAMS_FOUND AND (NOT WIN32 OR Boost_ZLIB_FOUND))
    qt5_wrap_ui( classificationUI_FILES Classification_widget.ui Classification_advanced_widget.ui )
    polyhedron_demo_plugin(classification_plugin Classification_plugin Point_set_item_classification Cluster_classification Surface_mesh_item_classification ${classificationUI_FILES})

    if(TARGET Boost::serialization AND TARGET Boost::iostreams)
      set(classification_linked_libraries ${classification_linked_libraries}
        Boost::serialization Boost::iostreams)
    else()
      set(classification_linked_libraries ${classification_linked_libraries}
        ${Boost_SERIALIZATION_LIBRARY}
        ${Boost_IOSTREAMS_LIBRARY})
    endif()

    find_package(OpenCV QUIET COMPONENTS core ml) # Need core + machine learning
    if (OpenCV_FOUND)
      message(STATUS "Found OpenCV ${OpenCV_VERSION}")
      include_directories(${OpenCV_INCLUDE_DIRS})
      set(classification_linked_libraries ${classification_linked_libraries}
        ${OpenCV_LIBS})
      set(classification_compile_definitions ${classification_compile_definitions}
        "-DCGAL_LINKED_WITH_OPENCV")
    else()
      message(STATUS "NOTICE: OpenCV was not found. OpenCV random forest predicate for classification won't be available.")
    endif()
    
    target_link_libraries(classification_plugin PUBLIC scene_points_with_normal_item
      scene_polylines_item scene_polygon_soup_item scene_surface_mesh_item scene_surface_mesh_selection_item scene_color_ramp ${classification_linked_libraries})
    target_compile_definitions(classification_plugin PUBLIC ${classification_compile_definitions})
  else()
    message(STATUS "NOTICE: Boost Serialization or IO Streams or ZLIB not found. Classification plugin won't be available.")
  endif()


else(EIGEN3_FOUND)
  message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Classification plugin will not be available.")
endif()


