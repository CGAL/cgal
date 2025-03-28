include(CGALlab_macros)

if(TARGET CGAL::Eigen3_support)
  cgal_lab_plugin(jet_fitting_plugin Jet_fitting_plugin)
  target_link_libraries(
    jet_fitting_plugin PUBLIC scene_surface_mesh_item scene_polylines_item
                              CGAL::Eigen3_support)
else()
  message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Jet fitting plugin will not be available.")
endif()

if(TARGET CGAL::Eigen3_support)

  cgal_lab_plugin(interpolated_corrected_principal_curvatures_plugin Interpolated_corrected_principal_curvatures_plugin)
  target_link_libraries(
    interpolated_corrected_principal_curvatures_plugin PUBLIC scene_surface_mesh_item scene_polylines_item
                              CGAL::Eigen3_support)

else()
  message(
    STATUS
      "NOTICE: Eigen 3.1 (or greater) was not found. Interpolated corrected principal curvatures plugin will not be available."
  )
endif()

cgal_lab_plugin(extrude_plugin Extrude_plugin KEYWORDS PMP)
target_link_libraries(extrude_plugin PUBLIC scene_surface_mesh_item
                                            scene_selection_item)

if(TARGET CGAL::Eigen3_support)
  if("${Eigen3_VERSION}" VERSION_GREATER "3.1.90")
    qt6_wrap_ui( hole_fillingUI_FILES Hole_filling_widget.ui)
    cgal_lab_plugin(hole_filling_plugin Hole_filling_plugin ${hole_fillingUI_FILES} KEYWORDS PMP)
    target_link_libraries(hole_filling_plugin PUBLIC scene_surface_mesh_item scene_polylines_item scene_selection_item CGAL::Eigen3_support)

    qt6_wrap_ui( fairingUI_FILES Fairing_widget.ui)
    cgal_lab_plugin(fairing_plugin Fairing_plugin ${fairingUI_FILES} KEYWORDS PMP)
    target_link_libraries(fairing_plugin PUBLIC scene_selection_item CGAL::Eigen3_support)

    qt6_wrap_ui( Mean_curvature_flow_skeleton_pluginUI_FILES  Mean_curvature_flow_skeleton_plugin.ui)
    cgal_lab_plugin(mean_curvature_flow_skeleton_plugin Mean_curvature_flow_skeleton_plugin ${Mean_curvature_flow_skeleton_pluginUI_FILES})
    target_link_libraries(mean_curvature_flow_skeleton_plugin
      PUBLIC
        scene_surface_mesh_item
        scene_points_with_normal_item
        scene_polylines_item
        scene_mcf_item
        demo_framework
        CGAL::Eigen3_support)
    # The smoothing plugin can still do some things, even if Ceres is not found
    qt6_wrap_ui( smoothingUI_FILES Smoothing_plugin.ui Smoothing_tangential_relaxation.ui)
    cgal_lab_plugin(smoothing_plugin Smoothing_plugin ${smoothingUI_FILES})
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

    qt6_wrap_ui(remeshPlanarPatchesUI_FILES Remesh_planar_patches_dialog.ui)
    cgal_lab_plugin(remesh_planar_patches_plugin Remesh_planar_patches_plugin
                           ${remeshPlanarPatchesUI_FILES} KEYWORDS PMP)
    target_link_libraries(remesh_planar_patches_plugin PUBLIC scene_surface_mesh_item CGAL::Eigen3_support)

    if(CGAL_ENABLE_TESTING AND NOT CMAKE_VS_MSBUILD_COMMAND)
      set_tests_properties(
        "compilation of  extrude_plugin"
        "compilation of  fairing_plugin"
        "compilation of  hole_filling_plugin"
        "compilation of  remesh_planar_patches_plugin"
        "compilation of  smoothing_plugin"
      PROPERTIES RESOURCE_LOCK Selection_test_resources)
    endif()
  else()
    message(STATUS "NOTICE: The hole filling and fairing plugins require Eigen 3.2 (or higher) and will not be available.")
  endif()
else()
  message(STATUS "NOTICE: The hole filling and fairing plugins require Eigen 3.2 (or higher) and will not be available.")
endif()

qt6_wrap_ui(soupUI_FILES Repair_soup.ui)
cgal_lab_plugin(orient_soup_plugin Orient_soup_plugin ${soupUI_FILES}
                       KEYWORDS Classification PMP)
target_link_libraries(
  orient_soup_plugin PUBLIC scene_polygon_soup_item scene_surface_mesh_item
                            scene_polylines_item scene_points_with_normal_item)

cgal_lab_plugin(inside_out_plugin Inside_out_plugin KEYWORDS PMP)

target_link_libraries(inside_out_plugin PUBLIC scene_surface_mesh_item scene_polygon_soup_item)

cgal_lab_plugin(join_and_split_plugin Join_and_split_polyhedra_plugin KEYWORDS PMP)
target_link_libraries(join_and_split_plugin PUBLIC scene_surface_mesh_item scene_selection_item)

qt6_wrap_ui( point_inside_polyhedronUI_FILES Point_inside_polyhedron_widget.ui)
cgal_lab_plugin(point_inside_polyhedron_plugin Point_inside_polyhedron_plugin ${point_inside_polyhedronUI_FILES})
target_link_libraries(point_inside_polyhedron_plugin PUBLIC scene_surface_mesh_item scene_points_with_normal_item)

qt6_wrap_ui( polyhedron_slicerUI_FILES Polyhedron_slicer_widget.ui)
cgal_lab_plugin(polyhedron_slicer_plugin Polyhedron_slicer_plugin ${polyhedron_slicerUI_FILES})
target_link_libraries(polyhedron_slicer_plugin PUBLIC scene_surface_mesh_item scene_basic_objects scene_polylines_item)

cgal_lab_plugin(polyhedron_stitching_plugin Polyhedron_stitching_plugin KEYWORDS PMP)
target_link_libraries(polyhedron_stitching_plugin PUBLIC scene_surface_mesh_item scene_polylines_item)

qt6_wrap_ui( selectionUI_FILES Selection_widget.ui)
cgal_lab_plugin(selection_plugin Selection_plugin ${selectionUI_FILES} KEYWORDS PMP Viewer Classification Mesh_3)
target_link_libraries(selection_plugin PUBLIC scene_selection_item scene_points_with_normal_item scene_polylines_item)

#to keep it simple to compile
add_custom_target(self_intersection_plugin)
add_dependencies(self_intersection_plugin selection_plugin)

cgal_lab_plugin(triangulate_facets_plugin Triangulate_facets_plugin KEYWORDS PMP)
target_link_libraries(triangulate_facets_plugin PUBLIC scene_surface_mesh_item scene_selection_item scene_polygon_soup_item)

cgal_lab_plugin(corefinement_plugin Corefinement_plugin KEYWORDS PMP)
target_link_libraries(corefinement_plugin PUBLIC scene_surface_mesh_item)

cgal_lab_plugin(surface_intersection_plugin Surface_intersection_plugin
                       KEYWORDS PMP)
target_link_libraries(
  surface_intersection_plugin
  PUBLIC scene_surface_mesh_item scene_polylines_item
         scene_points_with_normal_item)

qt6_wrap_ui( repairUI_FILES RemoveNeedlesDialog.ui SelfSnapDialog.ui AddBboxDialog.ui)
cgal_lab_plugin(repair_polyhedron_plugin Repair_polyhedron_plugin  ${repairUI_FILES} KEYWORDS PMP)
target_link_libraries(repair_polyhedron_plugin PUBLIC scene_points_with_normal_item scene_surface_mesh_item scene_polygon_soup_item)
if(TARGET CGAL::TBB_support)
  target_link_libraries(repair_polyhedron_plugin PUBLIC CGAL::TBB_support)
endif()

if(TARGET CGAL::Eigen3_support)
  qt6_wrap_ui(isotropicRemeshingUI_FILES Isotropic_remeshing_dialog.ui)
  cgal_lab_plugin(isotropic_remeshing_plugin Isotropic_remeshing_plugin
          ${isotropicRemeshingUI_FILES} KEYWORDS PMP)
  target_link_libraries(isotropic_remeshing_plugin PUBLIC scene_surface_mesh_item
          scene_selection_item CGAL::Eigen3_support)

  if(TARGET CGAL::TBB_support)
    target_link_libraries(isotropic_remeshing_plugin PUBLIC CGAL::TBB_support)
  endif()

else()
  message(STATUS "NOTICE: Eigen 3.1 (or greater) was not found. Isotropic remeshing plugin will not be available.")
endif()

cgal_lab_plugin(distance_plugin Distance_plugin KEYWORDS PMP)
target_link_libraries(distance_plugin PUBLIC scene_surface_mesh_item
                                             scene_color_ramp)
if(TARGET CGAL::TBB_support)
  target_link_libraries(distance_plugin PUBLIC CGAL::TBB_support)
endif()

cgal_lab_plugin(detect_sharp_edges_plugin Detect_sharp_edges_plugin KEYWORDS Viewer Mesh_3 PMP)
target_link_libraries(detect_sharp_edges_plugin PUBLIC scene_surface_mesh_item)

qt6_wrap_ui(randomPerturbationUI_FILES Random_perturbation_dialog.ui)
cgal_lab_plugin(random_perturbation_plugin Random_perturbation_plugin
                       ${randomPerturbationUI_FILES} KEYWORDS PMP)
target_link_libraries(random_perturbation_plugin PUBLIC scene_surface_mesh_item
                                                        scene_selection_item)

cgal_lab_plugin(degenerated_faces_plugin Degenerated_faces_plugin
                       KEYWORDS PMP)
target_link_libraries(degenerated_faces_plugin PUBLIC scene_surface_mesh_item
                                                      scene_selection_item)

qt6_wrap_ui(engravUI_FILES Engrave_dock_widget.ui)
cgal_lab_plugin(engrave_text_plugin Engrave_text_plugin
                       ${engravUI_FILES})
target_link_libraries(
  engrave_text_plugin PUBLIC scene_surface_mesh_item scene_selection_item
                             scene_polylines_item)

if(CGAL_ENABLE_TESTING AND NOT CMAKE_VS_MSBUILD_COMMAND)
  set_tests_properties(
    "compilation of  degenerated_faces_plugin"
    "compilation of  engrave_text_plugin"
    "compilation of  isotropic_remeshing_plugin"
    "compilation of  join_and_split_plugin"
    "compilation of  random_perturbation_plugin"
    "compilation of  selection_plugin"
    "compilation of  triangulate_facets_plugin"
    PROPERTIES RESOURCE_LOCK Selection_test_resources)
endif()
