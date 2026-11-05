include(CGALlab_macros)

cgal_lab_plugin(convex_hull_plugin Convex_hull_plugin KEYWORDS
                       PointSetProcessing)
target_link_libraries(
  convex_hull_plugin PRIVATE scene_points_with_normal_item scene_polylines_item
                             scene_selection_item scene_surface_mesh_item)


if(TEST "${compilation_of__convex_hull_plugin}")
  set_tests_properties(
    "${compilation_of__convex_hull_plugin}"
    PROPERTIES RESOURCE_LOCK Selection_test_resources)
endif()
