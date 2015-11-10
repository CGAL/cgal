include( polyhedron_demo_macros )
if ( EIGEN3_FOUND AND "${EIGEN3_VERSION}" VERSION_GREATER "3.1.90" )
  qt5_wrap_ui( editionUI_FILES Deform_mesh.ui )
  add_item(scene_edit_polyhedron_item Scene_edit_polyhedron_item.cpp
           ${editionUI_FILES})

  target_link_libraries(scene_edit_polyhedron_item scene_polyhedron_item scene_polyhedron_item_k_ring_selection)
  polyhedron_demo_plugin(edit_polyhedron_plugin Edit_polyhedron_plugin ${editionUI_FILES})
  target_link_libraries(edit_polyhedron_plugin scene_polyhedron_item scene_edit_polyhedron_item)
else()
  message(STATUS "NOTICE: The polyhedron edit plugin require Eigen 3.2 (or higher) and will not be available.")
endif()
