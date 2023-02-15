include(polyhedron_demo_macros)

remove_definitions(-DQT_STATICPLUGIN)

qt5_wrap_cpp(VOLUME_MOC_OUTFILES
             ${CMAKE_CURRENT_SOURCE_DIR}/Volume_plane_thread.h)
qt5_wrap_cpp(VOLUME_MOC_OUTFILES
             ${CMAKE_CURRENT_SOURCE_DIR}/Volume_plane_interface.h)

qt5_wrap_ui(tetRemeshingUI_FILES Tetrahedral_remeshing_dialog.ui)
polyhedron_demo_plugin(
  tetrahedral_remeshing_plugin Tetrahedral_remeshing_plugin
  ${tetRemeshingUI_FILES} KEYWORDS Tetrahedral_remeshing)
target_link_libraries(tetrahedral_remeshing_plugin PUBLIC scene_c3t3_item
                                                          ${OPENGL_gl_LIBRARY})
