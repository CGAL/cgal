include(polyhedron_demo_macros)

polyhedron_demo_plugin(gocad_plugin GOCAD_io_plugin KEYWORDS Viewer)
target_link_libraries(gocad_plugin PUBLIC scene_surface_mesh_item)

qt5_wrap_ui( funcUI_FILES Function_dialog.ui )
polyhedron_demo_plugin(io_implicit_function_plugin
  Implicit_function_io_plugin
  ${funcUI_FILES}
  KEYWORDS Viewer Mesh_3)
target_link_libraries(io_implicit_function_plugin PUBLIC scene_implicit_function_item)

polyhedron_demo_plugin(nef_io_plugin Nef_io_plugin KEYWORDS Viewer)
target_link_libraries(nef_io_plugin PUBLIC scene_nef_polyhedron_item)

polyhedron_demo_plugin(off_plugin OFF_io_plugin
  KEYWORDS Viewer Mesh_3 PointSetProcessing Classification PMP)
target_link_libraries(off_plugin PUBLIC scene_polygon_soup_item scene_points_with_normal_item scene_surface_mesh_item)

polyhedron_demo_plugin(off_to_nef_plugin OFF_to_nef_io_plugin KEYWORDS Viewer)
target_link_libraries(off_to_nef_plugin PUBLIC scene_nef_polyhedron_item)

polyhedron_demo_plugin(polylines_io_plugin Polylines_io_plugin KEYWORDS Viewer Mesh_3)
target_link_libraries(polylines_io_plugin PUBLIC scene_polylines_item)

polyhedron_demo_plugin(wkt_plugin WKT_io_plugin KEYWORDS Viewer PointSetProcessing Mesh_3)
target_link_libraries(wkt_plugin PUBLIC scene_polylines_item)

polyhedron_demo_plugin(stl_plugin STL_io_plugin KEYWORDS Viewer PMP)
target_link_libraries(stl_plugin PUBLIC scene_surface_mesh_item scene_polygon_soup_item)

polyhedron_demo_plugin(surf_io_plugin Surf_io_plugin KEYWORDS Viewer PMP)
target_link_libraries(surf_io_plugin PUBLIC scene_surface_mesh_item)

polyhedron_demo_plugin(lcc_io_plugin lcc_io_plugin KEYWORDS Viewer)
target_link_libraries(lcc_io_plugin PUBLIC scene_lcc_item)

find_package(VTK QUIET COMPONENTS vtkCommonCore vtkIOCore vtkIOLegacy vtkIOXML
                                  vtkFiltersCore vtkFiltersSources)
set_package_properties(
  VTK PROPERTIES
  DESCRIPTION "A library for image processing."
  PURPOSE "Can be used for I/O (DICOM, VTU, VTP.")

if(VTK_FOUND)
  if(VTK_USE_FILE)
    include(${VTK_USE_FILE})
  endif()
  if("${VTK_VERSION_MAJOR}" GREATER "5" OR VTK_VERSION VERSION_GREATER 5)
    if(TARGET VTK::CommonCore)
      set(VTK_LIBRARIES VTK::CommonCore VTK::IOCore VTK::IOLegacy VTK::IOXML
                        VTK::FiltersCore VTK::FiltersSources)
    endif()
    if(VTK_LIBRARIES)
      polyhedron_demo_plugin(vtk_plugin VTK_io_plugin KEYWORDS Viewer Mesh_3)
      target_link_libraries(vtk_plugin PUBLIC scene_surface_mesh_item scene_polylines_item scene_c3t3_item scene_points_with_normal_item
                                       ${VTK_LIBRARIES})

    else()
      message(
        STATUS
          "NOTICE : the vtk IO plugin needs VTK libraries and will not be compiled."
      )
    endif()
  else()
    message(
      STATUS
        "NOTICE : the vtk IO plugin needs VTK 6.0 or greater and will not be compiled (incorrect version found)."
    )
  endif()
else()
  message(
    STATUS
      "NOTICE : the vtk IO plugin needs VTK 6.0 or greater and will not be compiled."
  )
endif()
polyhedron_demo_plugin(xyz_plugin XYZ_io_plugin KEYWORDS Viewer PointSetProcessing Classification)
target_link_libraries(xyz_plugin PUBLIC scene_points_with_normal_item)

list(FIND CMAKE_CXX_COMPILE_FEATURES cxx_rvalue_references has_cxx_rvalues)
list(FIND CMAKE_CXX_COMPILE_FEATURES cxx_variadic_templates has_cxx_variadic)

if(has_cxx_rvalues LESS 0 OR has_cxx_variadic LESS 0)
  message(
    STATUS
      "NOTICE : LAS/PLY IO examples require a C++11 compiler and will not be compiled."
  )
else()
  set(needed_cxx_features cxx_rvalue_references cxx_variadic_templates)

  polyhedron_demo_plugin(ply_plugin PLY_io_plugin KEYWORDS Viewer PointSetProcessing Classification PMP)
  target_link_libraries(ply_plugin PUBLIC scene_points_with_normal_item scene_polygon_soup_item scene_surface_mesh_item scene_textured_item)
  target_compile_features(ply_plugin PRIVATE ${needed_cxx_features})

  if (TARGET CGAL::LASLIB_support)
    polyhedron_demo_plugin(las_plugin LAS_io_plugin KEYWORDS Viewer PointSetProcessing Classification)
    target_link_libraries(las_plugin PUBLIC scene_points_with_normal_item CGAL::LASLIB_support)
    target_compile_features(las_plugin PRIVATE ${needed_cxx_features})
  else()
    message(
      STATUS
        "NOTICE : the LAS IO plugin needs LAS libraries and will not be compiled."
    )
  endif()
endif()

find_path(
  3MF_INCLUDE_DIR
  NAMES Model/COM/NMR_DLLInterfaces.h
  PATH_SUFFIXES lib3mf
  DOC "Path to lib3MF headers")
find_library(
  3MF_LIBRARIES
  NAMES 3MF
  DOC "Path to the lib3MF library")

if(3MF_LIBRARIES
   AND 3MF_INCLUDE_DIR
   AND EXISTS "${3MF_INCLUDE_DIR}/Model/COM/NMR_DLLInterfaces.h")
  include_directories(${3MF_INCLUDE_DIR})
  polyhedron_demo_plugin(3mf_io_plugin 3mf_io_plugin KEYWORDS Viewer PMP)
  target_link_libraries(3mf_io_plugin PRIVATE scene_surface_mesh_item scene_points_with_normal_item scene_polylines_item ${3MF_LIBRARIES})
else()
  message(
    STATUS
      "NOTICE : The 3mf_io_plugin requires the lib3MF library in a version < 2.0, and will not be compiled."
  )
endif()
