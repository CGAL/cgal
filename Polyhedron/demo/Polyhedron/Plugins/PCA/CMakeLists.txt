include(polyhedron_demo_macros)
polyhedron_demo_plugin(pca_plugin Pca_plugin KEYWORDS PointSetProcessing)
target_link_libraries(
  pca_plugin PUBLIC scene_surface_mesh_item scene_points_with_normal_item
                    scene_basic_objects)

qt5_wrap_ui(transformUI_FILES Transformation_widget.ui MeshOnGrid_dialog.ui)

polyhedron_demo_plugin(affine_transform_plugin Affine_transform_plugin
                       ${transformUI_FILES} KEYWORDS PointSetProcessing)
target_link_libraries(
  affine_transform_plugin PUBLIC scene_surface_mesh_item scene_transform_item
                                 scene_points_with_normal_item)

polyhedron_demo_plugin(edit_box_plugin Edit_box_plugin)
target_link_libraries(edit_box_plugin PUBLIC scene_edit_box_item
                                             scene_surface_mesh_item)

qt5_wrap_ui(clipUI_FILES Clipping_box_widget.ui)
polyhedron_demo_plugin(clipping_box_plugin Clipping_box_plugin ${clipUI_FILES})
target_link_libraries(clipping_box_plugin PUBLIC scene_edit_box_item
                                                 scene_basic_objects)

polyhedron_demo_plugin(create_bbox_mesh_plugin Create_bbox_mesh_plugin)
target_link_libraries(create_bbox_mesh_plugin PUBLIC scene_surface_mesh_item)

polyhedron_demo_plugin(create_obb_mesh_plugin Create_obb_mesh_plugin)
target_link_libraries(
  create_obb_mesh_plugin PUBLIC scene_surface_mesh_item scene_selection_item
                                scene_points_with_normal_item)

qt5_wrap_ui(volumesUI_FILES Basic_generator_widget.ui)
polyhedron_demo_plugin(
  basic_generator_plugin Basic_generator_plugin ${volumesUI_FILES} KEYWORDS
  PolygonMesh PointSetProcessing)
target_link_libraries(
  basic_generator_plugin
  PUBLIC scene_surface_mesh_item scene_points_with_normal_item
         scene_polylines_item)
