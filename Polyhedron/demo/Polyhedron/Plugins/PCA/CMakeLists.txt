include( polyhedron_demo_macros )

if(${USE_PCH})
  add_library(PCA_headers pca_headers.cpp)
  target_precompile_headers(PCA_headers PUBLIC
  [["Scene_surface_mesh_item.h"]]
  [["Scene_points_with_normal_item.h"]]
  <QApplication>
  <CGAL/Three/Three.h>
  <QMainWindow>
  <QAction>
  <CGAL/Three/Viewer_interface.h>
  )

  if(TARGET CGAL::TBB_support)
    target_link_libraries(PCA_headers PUBLIC CGAL::TBB_support)
  endif()
  if(TARGET CGAL::Eigen_support)
    target_link_libraries(PCA_headers PUBLIC CGAL::Eigen_support)
  endif()

  if (TARGET CGAL::LASLIB_support)
    target_link_libraries(PCA_headers PUBLIC CGAL::LASLIB_support)
    if (MSVC)
      target_compile_definitions( PCA_headers PUBLIC "-D_CRT_SECURE_NO_DEPRECATE -D_CRT_SECURE_NO_WARNINGS")
    endif()
  endif()
endif()#use pch



polyhedron_demo_plugin(pca_plugin Pca_plugin KEYWORDS PointSetProcessing)
target_link_libraries(pca_plugin PUBLIC scene_surface_mesh_item scene_points_with_normal_item scene_basic_objects)
qt5_wrap_ui( transformUI_FILES Transformation_widget.ui MeshOnGrid_dialog.ui)
polyhedron_demo_plugin(affine_transform_plugin Affine_transform_plugin ${transformUI_FILES} KEYWORDS PMP PointSetProcessing)
target_link_libraries(affine_transform_plugin PUBLIC scene_surface_mesh_item scene_transform_item scene_points_with_normal_item)

polyhedron_demo_plugin(create_obb_mesh_plugin Create_obb_mesh_plugin)
target_link_libraries(create_obb_mesh_plugin PUBLIC scene_surface_mesh_item scene_selection_item scene_points_with_normal_item)

qt5_wrap_ui( volumesUI_FILES Basic_generator_widget.ui)
polyhedron_demo_plugin(basic_generator_plugin Basic_generator_plugin ${volumesUI_FILES} KEYWORDS PolygonMesh PointSetProcessing)
target_link_libraries(basic_generator_plugin PUBLIC scene_surface_mesh_item scene_points_with_normal_item scene_polylines_item)

if(${USE_PCH})
  foreach(tgt
  pca_plugin
  affine_transform_plugin
  create_obb_mesh_plugin
  basic_generator_plugin
  )
  target_precompile_headers(${tgt} REUSE_FROM PCA_headers)
  endforeach()
endif()#use pch


polyhedron_demo_plugin(create_bbox_mesh_plugin Create_bbox_mesh_plugin)
target_link_libraries(create_bbox_mesh_plugin PUBLIC scene_surface_mesh_item)
qt5_wrap_ui( clipUI_FILES Clipping_box_widget.ui )
polyhedron_demo_plugin(clipping_box_plugin Clipping_box_plugin ${clipUI_FILES})
target_link_libraries(clipping_box_plugin PUBLIC  scene_edit_box_item scene_basic_objects)
polyhedron_demo_plugin(edit_box_plugin Edit_box_plugin)
target_link_libraries(edit_box_plugin PUBLIC  scene_edit_box_item   scene_surface_mesh_item )
