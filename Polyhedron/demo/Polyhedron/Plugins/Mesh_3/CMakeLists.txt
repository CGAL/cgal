include( polyhedron_demo_macros )

if ( Boost_VERSION GREATER 103400 )
  remove_definitions(-DQT_STATICPLUGIN)

  qt5_wrap_cpp( VOLUME_MOC_OUTFILES ${CMAKE_CURRENT_SOURCE_DIR}/Volume_plane_thread.h )
  qt5_wrap_cpp( VOLUME_MOC_OUTFILES ${CMAKE_CURRENT_SOURCE_DIR}/Volume_plane_interface.h )
  qt5_wrap_ui( meshingUI_FILES  Meshing_dialog.ui Smoother_dialog.ui Local_optimizers_dialog.ui )
  qt5_generate_moc( "Polyhedron_demo_mesh_3_plugin_cgal_code.cpp" "${CMAKE_CURRENT_BINARY_DIR}/Scene_c3t3_item.moc" )
  polyhedron_demo_plugin(mesh_3_plugin Mesh_3_plugin
                         Mesh_3_plugin_cgal_code.cpp Meshing_thread.cpp split_polylines.cpp
                         ${meshingUI_FILES})
  target_link_libraries(mesh_3_plugin scene_polyhedron_item scene_polygon_soup_item scene_polylines_item scene_implicit_function_item scene_image_item
                        scene_surface_mesh_item scene_c3t3_item ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} )

  set(VTK_LIBS "")
  find_package(VTK QUIET COMPONENTS
    vtkImagingGeneral  vtkIOImage NO_MODULE)
  if (VTK_FOUND)
    include(${VTK_USE_FILE})
    if ("${VTK_VERSION_MAJOR}" GREATER "5")
      if(VTK_LIBRARIES)
        add_definitions(-DCGAL_USE_VTK)
        set(VTK_LIBS vtkImagingGeneral vtkIOImage)
      else()
        message(STATUS "NOTICE : the DICOM files (.dcm) need VTK libraries to be open and will not be able to.")
      endif()
    else()
      message(STATUS "NOTICE : the DICOM files (.dcm) need VTK libraries to be open and will not be able to.")
    endif()
  else()
    message(STATUS "NOTICE : the DICOM files (.dcm) need VTK libraries to be open and will not be able to.")
  endif()
  if(Boost_FILESYSTEM_FOUND)
    qt5_wrap_ui( imgUI_FILES Image_res_dialog.ui raw_image.ui)
    polyhedron_demo_plugin(io_image_plugin Io_image_plugin Volume_plane_intersection.cpp Raw_image_dialog.cpp ${imgUI_FILES} ${VOLUME_MOC_OUTFILES})
    target_link_libraries(io_image_plugin scene_image_item ${VTK_LIBS} )
  else()
      message( STATUS "NOTICE : the Io_image_plugin needs boost-filesystem to work and will not be compiled")
  endif()
  polyhedron_demo_plugin(mesh_3_optimization_plugin Optimization_plugin
    Optimization_plugin_cgal_code.cpp Optimizer_thread.cpp
    ${meshingUI_FILES})
  target_link_libraries(mesh_3_optimization_plugin scene_c3t3_item scene_polyhedron_item scene_surface_mesh_item scene_image_item scene_implicit_function_item )



  polyhedron_demo_plugin(c3t3_io_plugin C3t3_io_plugin)
  target_link_libraries(c3t3_io_plugin  scene_c3t3_item)
  else( Boost_VERSION GREATER 103400 )
    message(STATUS "warning: the plugin mesh_3_plugin requires Boost>=1.34.1 and will not be compiled.")
  endif( Boost_VERSION GREATER 103400 )
  qt5_wrap_ui( ribUI_FILES  Rib_dialog.ui)
  polyhedron_demo_plugin(c3t3_rib_exporter_plugin C3t3_rib_exporter_plugin ${ribUI_FILES})
  target_link_libraries(c3t3_rib_exporter_plugin scene_c3t3_item)
