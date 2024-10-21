include(CGALlab_macros)

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()
if(NOT CGAL_DISABLE_GMP)
  if(TARGET CGAL::Eigen3_support)
    find_package(CGAL REQUIRED COMPONENTS Core)

    qt6_wrap_ui(parameterizationUI_FILES Parameterization_widget.ui OTE_dialog.ui)
    cgal_lab_plugin(parameterization_plugin Parameterization_plugin
                           ${parameterizationUI_FILES})
    target_link_libraries(
      parameterization_plugin PRIVATE scene_surface_mesh_item scene_textured_item
                                      scene_selection_item CGAL::Eigen3_support)

    if(CGAL_ENABLE_TESTING AND NOT CMAKE_VS_MSBUILD_COMMAND)
      set_tests_properties(
        "compilation of  parameterization_plugin"
        PROPERTIES RESOURCE_LOCK Selection_test_resources)
    endif()

  else()
    message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. The Parameterization plugin will not be available.")
  endif()

  qt6_wrap_ui(segmentationUI_FILES Mesh_segmentation_widget.ui)
  cgal_lab_plugin(mesh_segmentation_plugin Mesh_segmentation_plugin
                         ${segmentationUI_FILES})
  target_link_libraries(mesh_segmentation_plugin PRIVATE scene_surface_mesh_item)

  if(CGAL_ENABLE_TESTING AND NOT CMAKE_VS_MSBUILD_COMMAND)
    set_tests_properties(
      "compilation of  mesh_segmentation_plugin"
      PROPERTIES RESOURCE_LOCK Selection_test_resources)
  endif()

  qt6_wrap_ui( mesh_simplificationUI_FILES  Mesh_simplification_dialog.ui)
  cgal_lab_plugin(mesh_simplification_plugin Mesh_simplification_plugin ${mesh_simplificationUI_FILES})
  target_link_libraries(mesh_simplification_plugin PRIVATE scene_surface_mesh_item scene_selection_item)

  qt6_wrap_ui(shortestPathUI_FILES Shortest_path_widget.ui)
  cgal_lab_plugin(shortest_path_plugin Shortest_path_plugin
                         ${shortestPathUI_FILES})
  target_link_libraries(shortest_path_plugin
      PRIVATE
      scene_surface_mesh_item
      scene_selection_item
      scene_shortest_path_item
      scene_basic_objects)

  qt6_wrap_ui(basicUI_FILES Surface_mesh_approximation_dockwidget.ui)
  cgal_lab_plugin(
    surface_mesh_approximation_plugin Surface_mesh_approximation_plugin
    ${basicUI_FILES} VSA_wrapper.cpp)
  target_link_libraries(
    surface_mesh_approximation_plugin
    PRIVATE scene_surface_mesh_item scene_polygon_soup_item scene_polylines_item CGAL::Eigen3_support)
else()
  message(STATUS "NOTICE: Some tests require the CGAL_Core library, and will not be compiled.")
endif()
