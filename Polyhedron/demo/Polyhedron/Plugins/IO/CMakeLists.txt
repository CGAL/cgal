include( polyhedron_demo_macros )

polyhedron_demo_plugin(gocad_plugin GOCAD_io_plugin)
target_link_libraries(gocad_plugin scene_polyhedron_item)

qt5_wrap_ui( funcUI_FILES Function_dialog.ui )
polyhedron_demo_plugin(io_implicit_function_plugin Implicit_function_io_plugin ${funcUI_FILES})
target_link_libraries(io_implicit_function_plugin scene_implicit_function_item)

polyhedron_demo_plugin(nef_io_plugin Nef_io_plugin)
target_link_libraries(nef_io_plugin  scene_nef_polyhedron_item)

qt5_wrap_ui( statisticsPolyhedronUI_FILES Polyhedron_demo_statistics_on_polyhedron_dialog.ui)
polyhedron_demo_plugin(off_plugin OFF_io_plugin)
target_link_libraries(off_plugin scene_polyhedron_item scene_polygon_soup_item scene_points_with_normal_item)

polyhedron_demo_plugin(off_to_nef_plugin OFF_to_nef_io_plugin)
target_link_libraries(off_to_nef_plugin scene_nef_polyhedron_item)

qt5_wrap_ui( polylines_ioUI_FILES Add_polylines_dialog.ui )
polyhedron_demo_plugin(polylines_io_plugin Polylines_io_plugin ${polylines_ioUI_FILES})
target_link_libraries(polylines_io_plugin scene_polylines_item)

polyhedron_demo_plugin(selection_io_plugin Selection_io_plugin)
target_link_libraries(selection_io_plugin scene_polyhedron_selection_item)

polyhedron_demo_plugin(stl_plugin STL_io_plugin)
target_link_libraries(stl_plugin scene_polyhedron_item scene_polygon_soup_item)

polyhedron_demo_plugin(surface_mesh_io_plugin Surface_mesh_io_plugin)
target_link_libraries(surface_mesh_io_plugin scene_surface_mesh_item scene_polygon_soup_item)



find_package(VTK QUIET COMPONENTS
  vtkCommonCore vtkIOCore vtkIOLegacy vtkIOXML vtkFiltersCore vtkFiltersSources)
if (VTK_FOUND)
  include(${VTK_USE_FILE})
  if ("${VTK_VERSION_MAJOR}" GREATER "5")
    if(VTK_LIBRARIES)
      polyhedron_demo_plugin(vtk_plugin VTK_io_plugin)
      target_link_libraries(vtk_plugin scene_polyhedron_item scene_polylines_item
                                        vtkCommonCore vtkIOCore vtkIOLegacy vtkIOXML
                                        vtkFiltersCore vtkFiltersSources)
    else()
      message(STATUS "NOTICE : the vtk IO plugin needs VTK libraries and will not be compiled.")
    endif()
  else()
    message(STATUS "NOTICE : the vtk IO plugin needs VTK 6.0 or greater and will not be compiled (incorrect version found).")
  endif()
else()
  message(STATUS "NOTICE : the vtk IO plugin needs VTK 6.0 or greater and will not be compiled.")
endif()
qt5_wrap_ui( xyzUI_FILES Add_point_set_dialog.ui)
polyhedron_demo_plugin(xyz_plugin XYZ_io_plugin ${xyzUI_FILES})
target_link_libraries(xyz_plugin scene_points_with_normal_item)

polyhedron_demo_plugin(ply_to_xyz_plugin PLY_to_xyz_io_plugin)
target_link_libraries(ply_to_xyz_plugin scene_points_with_normal_item)
