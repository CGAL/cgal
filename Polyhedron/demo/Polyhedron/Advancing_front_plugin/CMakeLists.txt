message(STATUS "Now CMAKE_MODULE_PATH is: ${CMAKE_MODULE_PATH}")
 include( polyhedron_demo_macros )
 
 qt5_wrap_ui( advancing_frontUI_FILES Polyhedron_demo_advancing_front_plugin.ui)
 polyhedron_demo_plugin(advancing_front_plugin Polyhedron_demo_advancing_front_plugin ${advancing_frontUI_FILES})
 target_link_libraries(advancing_front_plugin scene_polygon_soup_item scene_points_with_normal_item)

