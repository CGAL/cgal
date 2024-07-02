include(CGALlab_macros)

cgal_lab_plugin(convex_hull_plugin Convex_hull_plugin KEYWORDS
                       PointSetProcessing)
target_link_libraries(
  convex_hull_plugin PUBLIC scene_points_with_normal_item scene_polylines_item
                            scene_selection_item scene_surface_mesh_item)

cgal_lab_plugin(kernel_plugin Kernel_plugin)
target_link_libraries(kernel_plugin PUBLIC scene_surface_mesh_item)


if(CGAL_ENABLE_TESTING AND NOT CMAKE_VS_MSBUILD_COMMAND)
  set_tests_properties(
    "compilation of  convex_hull_plugin"
    PROPERTIES RESOURCE_LOCK Selection_test_resources)
endif()
