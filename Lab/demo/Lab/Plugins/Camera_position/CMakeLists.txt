
include( CGALlab_macros )
qt6_wrap_ui( cameraUI_FILES   Camera_positions_list.ui )
cgal_lab_plugin(camera_positions_plugin
    Camera_positions_plugin
    Camera_positions_list.cpp
    ${cameraUI_FILES}
    KEYWORDS Viewer)

add_dependencies(camera_positions_plugin demo_framework)
target_link_libraries(camera_positions_plugin PRIVATE demo_framework)
