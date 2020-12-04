include(polyhedron_demo_macros)
if(EIGEN3_FOUND AND "${EIGEN3_VERSION}" VERSION_GREATER "3.1.90")

  polyhedron_demo_plugin(edit_plugin Edit_polyhedron_plugin Deform_mesh.ui)
  target_link_libraries(edit_plugin PUBLIC scene_surface_mesh_item
                                           scene_edit_item scene_selection_item)
else()
  message(
    STATUS
      "NOTICE: The polyhedron edit plugin require Eigen 3.2 (or higher) and will not be available."
  )
endif()
