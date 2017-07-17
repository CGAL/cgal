include( polyhedron_demo_macros )

if(EIGEN3_FOUND)
  find_package(CGAL COMPONENTS Core)

  include(${CGAL_USE_FILE})

  if(CGAL_Core_FOUND)
    qt5_wrap_ui( parameterizationUI_FILES Parameterization_widget.ui ARAP_dialog.ui OTE_dialog.ui)
    polyhedron_demo_plugin(parameterization_plugin Parameterization_plugin ${parameterizationUI_FILES})
    target_link_libraries(parameterization_plugin ${CGAL_Core_LIBRARY} ${CGAL_Core_3RD_PARTY_LIBRARIES} scene_polyhedron_item scene_textured_polyhedron_item scene_polyhedron_selection_item)

    polyhedron_demo_plugin(parameterization_sm_plugin Parameterization_plugin ${parameterizationUI_FILES})
    target_link_libraries(parameterization_sm_plugin ${CGAL_Core_LIBRARY} ${CGAL_Core_3RD_PARTY_LIBRARIES} scene_surface_mesh_item scene_textured_surface_mesh_item scene_surface_mesh_selection_item)
    target_compile_definitions(parameterization_sm_plugin PUBLIC "-DUSE_SURFACE_MESH" )

  else(CGAL_Core_FOUND)
    message(STATUS "NOTICE: CGAL Core was not found. The Parameterization plugin will not be available.")
  endif(CGAL_Core_FOUND)

else(EIGEN3_FOUND)
  message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. The Parameterization plugin will not be available.")
endif(EIGEN3_FOUND)

qt5_wrap_ui( segmentationUI_FILES Mesh_segmentation_widget.ui)
polyhedron_demo_plugin(mesh_segmentation_plugin Mesh_segmentation_plugin ${segmentationUI_FILES})
target_link_libraries(mesh_segmentation_plugin scene_polyhedron_item scene_surface_mesh_item)

qt5_wrap_ui( mesh_simplificationUI_FILES  Mesh_simplification_dialog.ui)
polyhedron_demo_plugin(mesh_simplification_plugin Mesh_simplification_plugin ${mesh_simplificationUI_FILES})
target_link_libraries(mesh_simplification_plugin scene_polyhedron_item scene_polyhedron_selection_item)

polyhedron_demo_plugin(mesh_simplification_sm_plugin Mesh_simplification_plugin ${mesh_simplificationUI_FILES})
target_link_libraries(mesh_simplification_sm_plugin scene_surface_mesh_item scene_surface_mesh_selection_item)
target_compile_definitions(mesh_simplification_sm_plugin PUBLIC "-DUSE_SURFACE_MESH" )

qt5_wrap_ui( remeshingUI_FILES  Remeshing_dialog.ui)
polyhedron_demo_plugin(remeshing_plugin Remeshing_plugin Remeshing_plugin_cgal_code.cpp Mesher_base.cpp ${remeshingUI_FILES})
target_link_libraries(remeshing_plugin scene_polyhedron_item scene_polygon_soup_item scene_c2t3_item)

qt5_wrap_ui( shortestPathUI_FILES Shortest_path_widget.ui )
polyhedron_demo_plugin(shortest_path_plugin Shortest_path_plugin ${shortestPathUI_FILES})
target_link_libraries(shortest_path_plugin scene_polyhedron_item scene_polyhedron_shortest_path_item scene_basic_objects)

polyhedron_demo_plugin(shortest_path_sm_plugin Shortest_path_plugin ${shortestPathUI_FILES})
target_link_libraries(shortest_path_sm_plugin scene_surface_mesh_item scene_sm_shortest_path_item scene_basic_objects)
target_compile_definitions(shortest_path_sm_plugin PUBLIC "-DUSE_SURFACE_MESH" )
