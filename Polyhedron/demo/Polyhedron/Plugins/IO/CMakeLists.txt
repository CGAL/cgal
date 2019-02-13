include( polyhedron_demo_macros )

polyhedron_demo_plugin(gocad_plugin GOCAD_io_plugin KEYWORDS IO)
target_link_libraries(gocad_plugin PUBLIC scene_surface_mesh_item)

qt5_wrap_ui( funcUI_FILES Function_dialog.ui )
polyhedron_demo_plugin(io_implicit_function_plugin Implicit_function_io_plugin ${funcUI_FILES} KEYWORDS IO Mesh_3)
target_link_libraries(io_implicit_function_plugin PUBLIC scene_implicit_function_item)

polyhedron_demo_plugin(nef_io_plugin Nef_io_plugin KEYWORDS IO)
target_link_libraries(nef_io_plugin PUBLIC scene_nef_polyhedron_item)

polyhedron_demo_plugin(off_plugin OFF_io_plugin KEYWORDS IO Mesh_3 PointSetProcessing Classification)
target_link_libraries(off_plugin PUBLIC scene_polygon_soup_item scene_points_with_normal_item scene_surface_mesh_item)

polyhedron_demo_plugin(off_to_nef_plugin OFF_to_nef_io_plugin KEYWORDS IO)
target_link_libraries(off_to_nef_plugin PUBLIC scene_nef_polyhedron_item)

polyhedron_demo_plugin(polylines_io_plugin Polylines_io_plugin KEYWORDS IO Mesh_3)
target_link_libraries(polylines_io_plugin PUBLIC scene_polylines_item)

polyhedron_demo_plugin(stl_plugin STL_io_plugin KEYWORDS IO)
target_link_libraries(stl_plugin PUBLIC scene_surface_mesh_item scene_polygon_soup_item)


polyhedron_demo_plugin(surf_io_plugin Surf_io_plugin KEYWORDS IO)
target_link_libraries(surf_io_plugin PUBLIC scene_surface_mesh_item)



find_package(VTK QUIET COMPONENTS
  vtkCommonCore vtkIOCore vtkIOLegacy vtkIOXML vtkFiltersCore vtkFiltersSources)
if (VTK_FOUND)
  include(${VTK_USE_FILE})
  if ("${VTK_VERSION_MAJOR}" GREATER "5")
    if(VTK_LIBRARIES)
      polyhedron_demo_plugin(vtk_plugin VTK_io_plugin KEYWORDS IO Mesh_3)
      target_link_libraries(vtk_plugin PUBLIC scene_surface_mesh_item scene_polylines_item scene_c3t3_item scene_points_with_normal_item
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
polyhedron_demo_plugin(xyz_plugin XYZ_io_plugin KEYWORDS IO PointSetProcessing Classification)
target_link_libraries(xyz_plugin PUBLIC scene_points_with_normal_item)

list(FIND CMAKE_CXX_COMPILE_FEATURES cxx_rvalue_references has_cxx_rvalues)
list(FIND CMAKE_CXX_COMPILE_FEATURES cxx_variadic_templates has_cxx_variadic)

if(has_cxx_rvalues LESS 0 OR has_cxx_variadic LESS 0)
  message(STATUS "NOTICE: LAS/PLY IO examples require a C++11 compiler and will not be compiled.")
else()
  set(needed_cxx_features cxx_rvalue_references cxx_variadic_templates)
  
  polyhedron_demo_plugin(ply_plugin PLY_io_plugin KEYWORDS IO PointSetProcessing Classification)
  target_link_libraries(ply_plugin PUBLIC scene_points_with_normal_item scene_polygon_soup_item scene_surface_mesh_item scene_polygon_soup_item)
  target_compile_features(ply_plugin PRIVATE ${needed_cxx_features})
  
  if (LASLIB_FOUND)
    polyhedron_demo_plugin(las_plugin LAS_io_plugin KEYWORDS IO PointSetProcessing Classification)
    target_link_libraries(las_plugin PUBLIC scene_points_with_normal_item ${LASLIB_LIBRARIES})
    target_compile_features(las_plugin PRIVATE ${needed_cxx_features})
  else()
    message(STATUS "NOTICE : the LAS IO plugin needs LAS libraries and will not be compiled.")
  endif()
endif()
