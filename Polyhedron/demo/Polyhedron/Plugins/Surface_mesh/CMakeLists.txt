include(polyhedron_demo_macros)

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

if(TARGET CGAL::Eigen_support)
  find_package(CGAL REQUIRED COMPONENTS Core)

  include(${CGAL_USE_FILE})

  qt5_wrap_ui(parameterizationUI_FILES Parameterization_widget.ui OTE_dialog.ui)
  polyhedron_demo_plugin(parameterization_plugin Parameterization_plugin
                         ${parameterizationUI_FILES})
  target_link_libraries(
    parameterization_plugin PUBLIC scene_surface_mesh_item scene_textured_item
                                   scene_selection_item CGAL::Eigen_support)
else()
  message(
    STATUS
      "NOTICE: Eigen 3.1 (or greater) was not found. The Parameterization plugin will not be available."
  )
endif()

qt5_wrap_ui(segmentationUI_FILES Mesh_segmentation_widget.ui)
polyhedron_demo_plugin(mesh_segmentation_plugin Mesh_segmentation_plugin
                       ${segmentationUI_FILES})
target_link_libraries(mesh_segmentation_plugin PUBLIC scene_surface_mesh_item)


qt5_wrap_ui( mesh_simplificationUI_FILES  Mesh_simplification_dialog.ui)
polyhedron_demo_plugin(mesh_simplification_plugin Mesh_simplification_plugin ${mesh_simplificationUI_FILES})
target_link_libraries(mesh_simplification_plugin PUBLIC scene_surface_mesh_item scene_selection_item)

qt5_wrap_ui(remeshingUI_FILES Remeshing_dialog.ui)
polyhedron_demo_plugin(offset_meshing_plugin Offset_meshing_plugin
                       ${remeshingUI_FILES})
target_link_libraries(offset_meshing_plugin PUBLIC scene_surface_mesh_item
                                                   scene_polygon_soup_item)
if(TARGET CGAL::Eigen_support)
  target_link_libraries(offset_meshing_plugin PUBLIC CGAL::Eigen_support)
endif()
if(TARGET CGAL::TBB_support)
  target_link_libraries(offset_meshing_plugin PUBLIC CGAL::TBB_support)
endif()

qt5_wrap_ui(shortestPathUI_FILES Shortest_path_widget.ui)
polyhedron_demo_plugin(shortest_path_plugin Shortest_path_plugin
                       ${shortestPathUI_FILES})
target_link_libraries(
  shortest_path_plugin PUBLIC scene_surface_mesh_item scene_shortest_path_item
                              scene_basic_objects)

qt5_wrap_ui(basicUI_FILES Surface_mesh_approximation_dockwidget.ui)
polyhedron_demo_plugin(
  surface_mesh_approximation_plugin Surface_mesh_approximation_plugin
  ${basicUI_FILES} VSA_wrapper.cpp)
target_link_libraries(
  surface_mesh_approximation_plugin
  PUBLIC scene_surface_mesh_item scene_polygon_soup_item scene_polylines_item)
