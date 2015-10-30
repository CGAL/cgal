include( polyhedron_demo_macros )
qt5_wrap_ui( editionUI_FILES Deform_mesh.ui )
    #qt5_automoc( Scene_edit_polyhedron_item.cpp )
    add_item(scene_edit_polyhedron_item Scene_edit_polyhedron_item.cpp #Scene_edit_polyhedron_item.moc
      ${editionUI_FILES})

    target_link_libraries(scene_edit_polyhedron_item scene_polyhedron_item scene_polyhedron_item_k_ring_selection)

    polyhedron_demo_plugin(edit_polyhedron_plugin Polyhedron_demo_edit_polyhedron_plugin ${editionUI_FILES})
    target_link_libraries(edit_polyhedron_plugin scene_polyhedron_item scene_edit_polyhedron_item)