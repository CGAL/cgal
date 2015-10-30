include( polyhedron_demo_macros )
 if ( Boost_VERSION GREATER 103400 )
 qt5_wrap_ui( meshingUI_FILES  Meshing_dialog.ui )
    qt5_generate_moc( "Polyhedron_demo_mesh_3_plugin_cgal_code.cpp" "${CMAKE_CURRENT_BINARY_DIR}/Scene_c3t3_item.moc" )

    polyhedron_demo_plugin(mesh_3_plugin Polyhedron_demo_mesh_3_plugin 
      Polyhedron_demo_mesh_3_plugin_cgal_code.cpp #Scene_c3t3_item.moc
      ${meshingUI_FILES})
    target_link_libraries(mesh_3_plugin scene_polyhedron_item scene_polygon_soup_item ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})
  else( Boost_VERSION GREATER 103400 )
    message(STATUS "warning: the plugin mesh_3_plugin requires Boost>=1.34.1 and will not be compiled.")
  endif( Boost_VERSION GREATER 103400 )