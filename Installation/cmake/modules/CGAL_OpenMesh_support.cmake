
if(OpenMesh_FOUND AND NOT TARGET CGAL::OpenMesh_support)

  add_library(CGAL::OpenMesh_support INTERFACE IMPORTED)

  if(TARGET OpenMeshCore)
    target_link_libraries(CGAL::OpenMesh_support INTERFACE OpenMeshCore)
  endif()

  if(TARGET OpenMeshTools)
    target_link_libraries(CGAL::OpenMesh_support INTERFACE OpenMeshTools)
  endif()

  target_compile_definitions(CGAL::OpenMesh_support
    INTERFACE  "CGAL_USE_OPENMESH;NOMINMAX;_USE_MATH_DEFINES")

endif()
