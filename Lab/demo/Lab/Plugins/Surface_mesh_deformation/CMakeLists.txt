include(CGALlab_macros)

if(TARGET CGAL::Eigen3_support AND "${Eigen3_VERSION}" VERSION_GREATER "3.1.90")

  qt6_wrap_ui(editionUI_FILES Deform_mesh.ui)
  add_item(scene_edit_item
           Scene_edit_polyhedron_item.cpp
           ${editionUI_FILES})
  target_link_libraries(
    scene_edit_item PRIVATE CGAL::Eigen3_support scene_surface_mesh_item
                            scene_k_ring_selection scene_basic_objects)

  cgal_lab_plugin(edit_plugin Edit_polyhedron_plugin Deform_mesh.ui)
  target_link_libraries(edit_plugin PRIVATE scene_surface_mesh_item
                                            scene_edit_item scene_selection_item)

  if(CGAL_ENABLE_TESTING AND NOT CMAKE_VS_MSBUILD_COMMAND)
    set_tests_properties(
      "compilation of  edit_plugin"
      PROPERTIES RESOURCE_LOCK Selection_test_resources)
  endif()

else()
  message(STATUS "NOTICE: The polyhedron edit plugin requires Eigen 3.2 (or higher) and will not be available.")
endif()
