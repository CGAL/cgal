if(libpointmatcher_FOUND AND NOT TARGET CGAL::pointmatcher_support)
  add_library(CGAL::pointmatcher_support INTERFACE IMPORTED)
  set_target_properties(CGAL::pointmatcher_support PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "CGAL_LINKED_WITH_POINTMATCHER"
    INTERFACE_INCLUDE_DIRECTORIES "${libpointmatcher_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${libpointmatcher_LIBRARIES}")
endif()
