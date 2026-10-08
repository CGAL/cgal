include(CGALlab_macros)

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)


if(TARGET CGAL::Eigen3_support)
  cgal_lab_plugin(hexmeshing_lcc_plugin Hexmeshing_lcc_plugin KEYWORDS Viewer)
  target_link_libraries(hexmeshing_lcc_plugin PRIVATE CGAL::Eigen3_support scene_lcc_item scene_surface_mesh_item)
else()
  message(STATUS "NOTICE: The Hexmeshing_lcc_plugin requires Eigen3, which was not found. It will not be available.")
endif()

