include(CGALlab_macros)

if (NOT MSVC_VERSION OR MSVC_VERSION GREATER_EQUAL 1919 OR MSVC_VERSION LESS 1910)
  find_package(LASLIB QUIET)
  set_package_properties(
    LASLIB PROPERTIES
    DESCRIPTION "A library for LIDAR I/O."
    PURPOSE "Required for reading or writing LAS files.")
else()
  message(STATUS "NOTICE: the LAS reader does not work with your version of Visual Studio 2017.")
endif()

include(CGAL_LASLIB_support)

cgal_lab_plugin(gocad_plugin GOCAD_io_plugin KEYWORDS Viewer)
target_link_libraries(gocad_plugin PRIVATE scene_surface_mesh_item)

qt6_wrap_ui( funcUI_FILES Function_dialog.ui )
cgal_lab_plugin(io_implicit_function_plugin
  Implicit_function_io_plugin
  ${funcUI_FILES}
  KEYWORDS Viewer Mesh_3)
target_link_libraries(io_implicit_function_plugin PRIVATE scene_implicit_function_item)

cgal_lab_plugin(nef_io_plugin Nef_io_plugin KEYWORDS Viewer)
target_link_libraries(nef_io_plugin PRIVATE scene_nef_polyhedron_item)

cgal_lab_plugin(off_plugin OFF_io_plugin
  KEYWORDS Viewer Mesh_3 PointSetProcessing Classification PMP IO_surface_meshes)
target_link_libraries(off_plugin PRIVATE scene_polygon_soup_item scene_points_with_normal_item scene_surface_mesh_item)

cgal_lab_plugin(off_to_nef_plugin OFF_to_nef_io_plugin KEYWORDS Viewer)
target_link_libraries(off_to_nef_plugin PRIVATE scene_nef_polyhedron_item)

cgal_lab_plugin(polylines_io_plugin Polylines_io_plugin KEYWORDS Viewer Mesh_3)
target_link_libraries(polylines_io_plugin PRIVATE scene_polylines_item)

cgal_lab_plugin(wkt_plugin WKT_io_plugin KEYWORDS Viewer PointSetProcessing Mesh_3)
target_link_libraries(wkt_plugin PRIVATE scene_polylines_item)

cgal_lab_plugin(stl_plugin STL_io_plugin KEYWORDS Viewer PMP IO_surface_meshes)
target_link_libraries(stl_plugin PRIVATE scene_surface_mesh_item scene_polygon_soup_item)

cgal_lab_plugin(surf_io_plugin Surf_io_plugin KEYWORDS Viewer PMP IO_surface_meshes)
target_link_libraries(surf_io_plugin PRIVATE scene_surface_mesh_item)

cgal_lab_plugin(lcc_io_plugin lcc_io_plugin KEYWORDS Viewer)
target_link_libraries(lcc_io_plugin PRIVATE scene_lcc_item)

find_package(OpenMesh QUIET)
if(OpenMesh_FOUND)
  include(CGAL_OpenMesh_support)
  cgal_lab_plugin(om_plugin OM_io_plugin KEYWORDS Viewer PMP IO_surface_meshes)
  target_link_libraries(om_plugin PUBLIC scene_surface_mesh_item scene_polygon_soup_item scene_selection_item)
  target_link_libraries(om_plugin PRIVATE CGAL::OpenMesh_support)
else()
  message(STATUS "NOTICE: the OM IO plugin needs OpenMesh libraries and will not be compiled.")
endif()

set(CMAKE_POLICY_DEFAULT_CMP0167 NEW)
find_package(VTK 9.0 QUIET COMPONENTS CommonCore IOCore IOLegacy IOXML FiltersCore FiltersSources)
set_package_properties(
  VTK PROPERTIES
  DESCRIPTION "A library for image processing."
  PURPOSE "Can be used for I/O (DICOM, VTU, VTP).")

cgal_lab_plugin(triangulation_3_io_plugin triangulation_3_io_plugin KEYWORDS Viewer)
target_link_libraries(triangulation_3_io_plugin PRIVATE scene_triangulation_3_item)
if(TARGET CGAL::TBB_support)
  target_link_libraries(triangulation_3_io_plugin PRIVATE CGAL::TBB_support)
endif()

if(VTK_FOUND AND VTK_LIBRARIES)
  message(STATUS "VTK ${VTK_VERSION} found ${VTK_LIBRARIES}")
  cgal_lab_plugin(vtk_plugin VTK_io_plugin KEYWORDS Viewer IO_surface_meshes)
  target_link_libraries(vtk_plugin
      PRIVATE
      scene_surface_mesh_item
      scene_polygon_soup_item
      scene_polylines_item
      scene_c3t3_item
      scene_points_with_normal_item
      ${VTK_LIBRARIES}
  )
  target_compile_definitions(vtk_plugin  PRIVATE -DCGAL_USE_VTK -DNOMINMAX)
  CGAL_Lab_target_use_c3t3_type(vtk_plugin)
  set_target_properties(vtk_plugin PROPERTIES DISABLE_PRECOMPILE_HEADERS ON)
else()
  message(STATUS "NOTICE: the vtk IO plugin needs VTK 9.0 or greater and will not be compiled.")
endif()
cgal_lab_plugin(xyz_plugin XYZ_io_plugin KEYWORDS Viewer PointSetProcessing Classification)
target_link_libraries(xyz_plugin PRIVATE scene_points_with_normal_item)

cgal_lab_plugin(ply_plugin PLY_io_plugin KEYWORDS Viewer PointSetProcessing Classification PMP IO_surface_meshes)
target_link_libraries(ply_plugin PRIVATE scene_points_with_normal_item scene_polygon_soup_item scene_surface_mesh_item scene_textured_item)
if (TARGET CGAL::LASLIB_support)
  cgal_lab_plugin(las_plugin LAS_io_plugin KEYWORDS Viewer PointSetProcessing Classification)
  target_link_libraries(las_plugin PRIVATE scene_points_with_normal_item CGAL::LASLIB_support)
  if(MSVC)
    target_compile_definitions(
      las_plugin
      PRIVATE "-D_CRT_SECURE_NO_DEPRECATE -D_CRT_SECURE_NO_WARNINGS")
  endif()
else()
  message(STATUS "NOTICE: the LAS IO plugin needs LAS libraries and will not be compiled.")
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
  cgal_lab_plugin(io_3mf_plugin 3mf_io_plugin KEYWORDS Viewer PMP IO_surface_meshes)
  target_link_libraries(io_3mf_plugin PRIVATE scene_surface_mesh_item scene_points_with_normal_item scene_polylines_item ${3MF_LIBRARIES})
  target_compile_definitions(io_3mf_plugin PRIVATE -DCGAL_LINKED_WITH_3MF)
else()
  message(STATUS "NOTICE: The 3mf_io_plugin requires the lib3MF library in a version < 2.0, and will not be compiled.")
endif()
