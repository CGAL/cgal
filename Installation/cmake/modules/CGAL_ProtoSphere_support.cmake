if(ProtoSphere_FOUND AND NOT TARGET CGAL::ProtoSphere_support)
  add_library(CGAL::ProtoSphere_support INTERFACE IMPORTED)
  set_target_properties(CGAL::ProtoSphere_support PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${ProtoSphere_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "ProtoSphere::ProtoSphereExtended")
endif()
