include(CGALlab_macros)
cgal_lab_plugin(pca_plugin Pca_plugin KEYWORDS PointSetProcessing)
target_link_libraries(
  pca_plugin PRIVATE scene_surface_mesh_item scene_points_with_normal_item
                     scene_basic_objects CGAL::Eigen3_support)

qt6_wrap_ui(transformUI_FILES Transformation_widget.ui MeshOnGrid_dialog.ui)

cgal_lab_plugin(affine_transform_plugin Affine_transform_plugin
                       ${transformUI_FILES} KEYWORDS PointSetProcessing)
target_link_libraries(affine_transform_plugin PRIVATE scene_surface_mesh_item
                                                      scene_polygon_soup_item
                                                      scene_points_with_normal_item
                                                      scene_aff_transformed_item
                                                      scene_aff_transformed_point_set_item
                                                      scene_aff_transformed_polygon_soup_item
                                                      scene_aff_transformed_surface_mesh_item)

cgal_lab_plugin(edit_box_plugin Edit_box_plugin)
target_link_libraries(edit_box_plugin PRIVATE scene_edit_box_item
                                              scene_surface_mesh_item)

qt6_wrap_ui(clipUI_FILES Clipping_box_widget.ui)
cgal_lab_plugin(clipping_box_plugin Clipping_box_plugin ${clipUI_FILES})
target_link_libraries(clipping_box_plugin PRIVATE scene_edit_box_item
                                                  scene_basic_objects)

cgal_lab_plugin(create_bbox_mesh_plugin Create_bbox_mesh_plugin)
target_link_libraries(create_bbox_mesh_plugin PRIVATE scene_surface_mesh_item)

cgal_lab_plugin(create_obb_mesh_plugin Create_obb_mesh_plugin)
target_link_libraries(create_obb_mesh_plugin PRIVATE scene_surface_mesh_item
                                                     scene_polygon_soup_item
                                                     scene_selection_item
                                                     scene_points_with_normal_item CGAL::Eigen3_support)

qt6_wrap_ui(volumesUI_FILES Basic_generator_widget.ui)
cgal_lab_plugin(
  basic_generator_plugin Basic_generator_plugin ${volumesUI_FILES} KEYWORDS
  PolygonMesh PointSetProcessing)
target_link_libraries(
  basic_generator_plugin
  PRIVATE scene_surface_mesh_item scene_points_with_normal_item
          scene_polylines_item)

if(CGAL_ENABLE_TESTING AND NOT CMAKE_VS_MSBUILD_COMMAND)
  set_tests_properties(
    "compilation of  create_obb_mesh_plugin"
    PROPERTIES RESOURCE_LOCK Selection_test_resources)
endif()
