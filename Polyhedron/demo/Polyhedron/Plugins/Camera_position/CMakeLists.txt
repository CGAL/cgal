include( polyhedron_demo_macros )
  qt5_wrap_ui( cameraUI_FILES   Camera_positions_list.ui )
  polyhedron_demo_plugin(camera_positions_plugin 
    Camera_positions_plugin
    Camera_positions_list 
    ${cameraUI_FILES})

  target_link_libraries( camera_positions_plugin ${QGLVIEWER_LIBRARIES} demo_framework) 

