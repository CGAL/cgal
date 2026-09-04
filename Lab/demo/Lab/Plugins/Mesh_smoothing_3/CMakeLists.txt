include(CGALlab_macros)

remove_definitions(-DQT_STATICPLUGIN)

qt6_wrap_cpp(VOLUME_MOC_OUTFILES
             ${CMAKE_CURRENT_SOURCE_DIR}/Volume_plane_thread.h)
qt6_wrap_cpp(VOLUME_MOC_OUTFILES
             ${CMAKE_CURRENT_SOURCE_DIR}/Volume_plane_interface.h)

cgal_lab_plugin(
  mesh_smoothing_3_plugin Mesh_smoothing_3_plugin
  KEYWORDS Tetrahedral_remeshing)
target_link_libraries(mesh_smoothing_3_plugin PRIVATE scene_c3t3_item
                                                           ${OPENGL_gl_LIBRARY})
CGAL_Lab_target_use_c3t3_type(mesh_smoothing_3_plugin)
