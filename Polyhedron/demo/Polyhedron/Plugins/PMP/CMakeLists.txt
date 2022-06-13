include(polyhedron_demo_macros)
if(TARGET CGAL::Eigen3_support)

  polyhedron_demo_plugin(jet_fitting_plugin Jet_fitting_plugin)
  target_link_libraries(
    jet_fitting_plugin PUBLIC scene_surface_mesh_item scene_polylines_item
                              CGAL::Eigen3_support)

else()
  message(
    STATUS
      "NOTICE: Eigen 3.1 (or greater) was not found. Jet fitting plugin will not be available."
  )
endif()

polyhedron_demo_plugin(extrude_plugin Extrude_plugin KEYWORDS PMP)
target_link_libraries(extrude_plugin PUBLIC scene_surface_mesh_item
                                            scene_selection_item)

if(TARGET CGAL::Eigen3_support)
  if("${EIGEN3_VERSION}" VERSION_GREATER "3.1.90")
    qt5_wrap_ui( hole_fillingUI_FILES Hole_filling_widget.ui)
    polyhedron_demo_plugin(hole_filling_plugin Hole_filling_plugin ${hole_fillingUI_FILES} KEYWORDS PMP)
    target_link_libraries(hole_filling_plugin PUBLIC scene_surface_mesh_item scene_polylines_item scene_selection_item CGAL::Eigen3_support)

    qt5_wrap_ui( fairingUI_FILES Fairing_widget.ui)
    polyhedron_demo_plugin(fairing_plugin Fairing_plugin ${fairingUI_FILES} KEYWORDS PMP)
    target_link_libraries(fairing_plugin PUBLIC scene_selection_item CGAL::Eigen3_support)

    qt5_wrap_ui( Mean_curvature_flow_skeleton_pluginUI_FILES  Mean_curvature_flow_skeleton_plugin.ui)
    polyhedron_demo_plugin(mean_curvature_flow_skeleton_plugin Mean_curvature_flow_skeleton_plugin ${Mean_curvature_flow_skeleton_pluginUI_FILES})
    target_link_libraries(mean_curvature_flow_skeleton_plugin
      PUBLIC
        scene_surface_mesh_item
        scene_points_with_normal_item
        scene_polylines_item
        scene_mcf_item
        demo_framework
        CGAL::Eigen3_support)
      # The smoothing plugin can still do some things, even if Ceres is not found
      qt5_wrap_ui( smoothingUI_FILES Smoothing_plugin.ui Smoothing_tangential_relaxation.ui)
      polyhedron_demo_plugin(smoothing_plugin Smoothing_plugin ${smoothingUI_FILES})
      target_link_libraries(smoothing_plugin PUBLIC scene_surface_mesh_item scene_selection_item CGAL::Eigen3_support)
      find_package(Ceres QUIET)
      include(CGAL_Ceres_support)
      if(TARGET CGAL::Ceres_support)
        target_link_libraries(smoothing_plugin PUBLIC CGAL::Ceres_support)
      endif()

      set_package_properties(
        Ceres PROPERTIES
        DESCRIPTION "A large scale non-linear optimization library."
        PURPOSE "Can be used as a solver in the smoothing plugin.")
        target_link_libraries(extrude_plugin PUBLIC CGAL::Eigen3_support)

        if(BUILD_TESTING AND NOT CMAKE_VS_MSBUILD_COMMAND)
          set_tests_properties(
            compilation_of__extrude_plugin
            compilation_of__fairing_plugin
            PROPERTIES RESOURCE_LOCK Selection_test_resources)
          set_tests_properties(
          compilation_of__hole_filling_plugin
          compilation_of__smoothing_plugin
          PROPERTIES RESOURCE_LOCK Selection_test_resources)
        endif()
  else()
    message(
      STATUS
        "NOTICE: The hole filling and fairing plugins require Eigen 3.2 (or higher) and will not be available."
    )
  endif()

else(EIGEN3_FOUND)
  message(
    STATUS
      "NOTICE: The hole filling and fairing plugins require Eigen 3.2 (or higher) and will not be available."
  )
endif()

qt5_wrap_ui(soupUI_FILES Repair_soup.ui)
polyhedron_demo_plugin(orient_soup_plugin Orient_soup_plugin ${soupUI_FILES}
                       KEYWORDS Classification PMP)
target_link_libraries(
  orient_soup_plugin PUBLIC scene_polygon_soup_item scene_surface_mesh_item
                            scene_polylines_item scene_points_with_normal_item)

polyhedron_demo_plugin(inside_out_plugin Inside_out_plugin KEYWORDS PMP)

target_link_libraries(inside_out_plugin PUBLIC scene_surface_mesh_item scene_polygon_soup_item)

polyhedron_demo_plugin(join_and_split_plugin Join_and_split_polyhedra_plugin KEYWORDS PMP)
target_link_libraries(join_and_split_plugin PUBLIC scene_surface_mesh_item scene_selection_item)

qt5_wrap_ui( point_inside_polyhedronUI_FILES Point_inside_polyhedron_widget.ui)
polyhedron_demo_plugin(point_inside_polyhedron_plugin Point_inside_polyhedron_plugin ${point_inside_polyhedronUI_FILES})
target_link_libraries(point_inside_polyhedron_plugin PUBLIC scene_surface_mesh_item scene_points_with_normal_item)

qt5_wrap_ui( polyhedron_slicerUI_FILES Polyhedron_slicer_widget.ui)
polyhedron_demo_plugin(polyhedron_slicer_plugin Polyhedron_slicer_plugin ${polyhedron_slicerUI_FILES})
target_link_libraries(polyhedron_slicer_plugin PUBLIC scene_surface_mesh_item scene_basic_objects scene_polylines_item)

polyhedron_demo_plugin(polyhedron_stitching_plugin Polyhedron_stitching_plugin KEYWORDS PMP)
target_link_libraries(polyhedron_stitching_plugin PUBLIC scene_surface_mesh_item scene_polylines_item)

qt5_wrap_ui( selectionUI_FILES Selection_widget.ui)
polyhedron_demo_plugin(selection_plugin Selection_plugin ${selectionUI_FILES} KEYWORDS PMP Viewer Classification Mesh_3)
target_link_libraries(selection_plugin PUBLIC scene_selection_item scene_points_with_normal_item scene_polylines_item)

#to keep it simple to compile
add_custom_target(self_intersection_plugin)
add_dependencies(self_intersection_plugin selection_plugin)

polyhedron_demo_plugin(triangulate_facets_plugin Triangulate_facets_plugin KEYWORDS PMP)
target_link_libraries(triangulate_facets_plugin PUBLIC scene_surface_mesh_item scene_selection_item)

polyhedron_demo_plugin(corefinement_plugin Corefinement_plugin KEYWORDS PMP)
target_link_libraries(corefinement_plugin PUBLIC scene_surface_mesh_item)

polyhedron_demo_plugin(surface_intersection_plugin Surface_intersection_plugin
                       KEYWORDS PMP)
target_link_libraries(
  surface_intersection_plugin
  PUBLIC scene_surface_mesh_item scene_polylines_item
         scene_points_with_normal_item)

qt5_wrap_ui( repairUI_FILES RemoveNeedlesDialog.ui SelfSnapDialog.ui)
polyhedron_demo_plugin(repair_polyhedron_plugin Repair_polyhedron_plugin  ${repairUI_FILES} KEYWORDS PMP)
target_link_libraries(repair_polyhedron_plugin PUBLIC scene_points_with_normal_item scene_surface_mesh_item)

qt5_wrap_ui(isotropicRemeshingUI_FILES Isotropic_remeshing_dialog.ui)
polyhedron_demo_plugin(isotropic_remeshing_plugin Isotropic_remeshing_plugin
                       ${isotropicRemeshingUI_FILES} KEYWORDS PMP)
target_link_libraries(isotropic_remeshing_plugin PUBLIC scene_surface_mesh_item
                                                        scene_selection_item)

if(TARGET CGAL::TBB_support)
  target_link_libraries(isotropic_remeshing_plugin PUBLIC CGAL::TBB_support)
endif()

polyhedron_demo_plugin(distance_plugin Distance_plugin KEYWORDS PMP)
target_link_libraries(distance_plugin PUBLIC scene_surface_mesh_item
                                             scene_color_ramp)
if(TARGET CGAL::TBB_support)
  target_link_libraries(distance_plugin PUBLIC CGAL::TBB_support)
endif()

polyhedron_demo_plugin(detect_sharp_edges_plugin Detect_sharp_edges_plugin KEYWORDS Viewer Mesh_3 PMP)
target_link_libraries(detect_sharp_edges_plugin PUBLIC scene_surface_mesh_item)

qt5_wrap_ui(randomPerturbationUI_FILES Random_perturbation_dialog.ui)
polyhedron_demo_plugin(random_perturbation_plugin Random_perturbation_plugin
                       ${randomPerturbationUI_FILES} KEYWORDS PMP)
target_link_libraries(random_perturbation_plugin PUBLIC scene_surface_mesh_item
                                                        scene_selection_item)

polyhedron_demo_plugin(degenerated_faces_plugin Degenerated_faces_plugin
                       KEYWORDS PMP)
target_link_libraries(degenerated_faces_plugin PUBLIC scene_surface_mesh_item
                                                      scene_selection_item)

qt5_wrap_ui(engravUI_FILES Engrave_dock_widget.ui)
polyhedron_demo_plugin(engrave_text_plugin Engrave_text_plugin
                       ${engravUI_FILES})
target_link_libraries(
  engrave_text_plugin PUBLIC scene_surface_mesh_item scene_selection_item
                             scene_polylines_item)

if(BUILD_TESTING AND NOT CMAKE_VS_MSBUILD_COMMAND)
  set_tests_properties(
    compilation_of__join_and_split_plugin
    compilation_of__selection_plugin
    compilation_of__triangulate_facets_plugin
    compilation_of__isotropic_remeshing_plugin
    compilation_of__random_perturbation_plugin
    compilation_of__engrave_text_plugin
    compilation_of__degenerated_faces_plugin
    PROPERTIES RESOURCE_LOCK Selection_test_resources)
endif()
