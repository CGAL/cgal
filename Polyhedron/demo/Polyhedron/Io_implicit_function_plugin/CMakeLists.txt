include( polyhedron_demo_macros )
  qt5_wrap_ui( funcUI_FILES Function_dialog.ui )
  polyhedron_demo_plugin(io_implicit_function_plugin Io_implicit_function_plugin ${funcUI_FILES})
  target_link_libraries(io_implicit_function_plugin scene_implicit_function_item)