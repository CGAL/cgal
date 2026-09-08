include(CGALlab_macros)

if(TARGET CGAL::Eigen3_support)

  cgal_lab_plugin(
    mesh_smoothing_3_plugin Mesh_smoothing_3_plugin
    KEYWORDS Mesh_smoothing)
  target_link_libraries(mesh_smoothing_3_plugin PRIVATE scene_c3t3_item CGAL::Eigen3_support
                                                            ${OPENGL_gl_LIBRARY})
  CGAL_Lab_target_use_c3t3_type(mesh_smoothing_3_plugin)

  if(TARGET CGAL::TBB_support)
    target_link_libraries(mesh_smoothing_3_plugin PRIVATE CGAL::TBB_support)
  endif()

  find_package(OpenMP COMPONENTS CXX QUIET)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(mesh_smoothing_3_plugin PRIVATE OpenMP::OpenMP_CXX)
  endif()

else()
  message(STATUS "NOTICE: The mesh smoothing 3 plugin requires Eigen 3.3 (or higher) and will not be available.")
endif()
