include(polyhedron_demo_macros)
qt5_wrap_ui(clip_polyhedronUI_FILES Clip_polyhedron_plugin.ui)
polyhedron_demo_plugin(clip_polyhedron_plugin Clip_polyhedron_plugin
                       ${clip_polyhedronUI_FILES})
target_link_libraries(
  clip_polyhedron_plugin PUBLIC scene_surface_mesh_item scene_basic_objects
                                scene_movable_sm_item)

polyhedron_demo_plugin(point_set_from_vertices_plugin
                       Point_set_from_vertices_plugin)
target_link_libraries(
  point_set_from_vertices_plugin
  PUBLIC scene_surface_mesh_item scene_polygon_soup_item
         scene_points_with_normal_item)

polyhedron_demo_plugin(diff_between_meshes_plugin Diff_between_meshes_plugin)
target_link_libraries(diff_between_meshes_plugin PUBLIC scene_surface_mesh_item)

 qt5_wrap_ui( animateUI_FILES Animate_widget.ui )
 polyhedron_demo_plugin(animate_mesh_plugin Animate_mesh_plugin ${animateUI_FILES})
 target_link_libraries(animate_mesh_plugin PUBLIC scene_surface_mesh_item)

 if( METIS_FOUND )
   qt5_wrap_ui( partitionUI_FILES PartitionDialog.ui )
   polyhedron_demo_plugin(partition_plugin Partition_graph_plugin ${partitionUI_FILES})
   target_link_libraries(partition_plugin PUBLIC scene_surface_mesh_item ${METIS_LIBRARIES} )
else()
  message(
    "NOTICE : the Partition plugin needs METIS libraries and will not be compiled."
  )
endif()
