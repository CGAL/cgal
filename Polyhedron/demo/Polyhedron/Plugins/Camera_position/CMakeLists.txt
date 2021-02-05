
include( polyhedron_demo_macros )
qt5_wrap_ui( cameraUI_FILES   Camera_positions_list.ui )
polyhedron_demo_plugin(camera_positions_plugin
    Camera_positions_plugin
    Camera_positions_list.cpp
    ${cameraUI_FILES}
    KEYWORDS Viewer)

target_link_libraries( camera_positions_plugin PUBLIC demo_framework)
target_link_libraries(camera_positions_plugin PUBLIC demo_framework)
