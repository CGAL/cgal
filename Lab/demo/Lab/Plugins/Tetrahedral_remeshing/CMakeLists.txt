include(CGALlab_macros)

remove_definitions(-DQT_STATICPLUGIN)

qt6_wrap_cpp(VOLUME_MOC_OUTFILES
             ${CMAKE_CURRENT_SOURCE_DIR}/Volume_plane_thread.h)
qt6_wrap_cpp(VOLUME_MOC_OUTFILES
             ${CMAKE_CURRENT_SOURCE_DIR}/Volume_plane_interface.h)

qt6_wrap_ui(tetRemeshingUI_FILES Tetrahedral_remeshing_dialog.ui)
cgal_lab_plugin(
  tetrahedral_remeshing_plugin Tetrahedral_remeshing_plugin
  ${tetRemeshingUI_FILES} KEYWORDS Tetrahedral_remeshing)
target_link_libraries(tetrahedral_remeshing_plugin PRIVATE scene_c3t3_item
                                                           ${OPENGL_gl_LIBRARY})
CGAL_Lab_target_use_c3t3_type(tetrahedral_remeshing_plugin)
find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  target_link_libraries(tetrahedral_remeshing_plugin PRIVATE CGAL::Eigen3_support)
endif()
