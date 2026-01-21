if(NanoSVG_FOUND AND NOT TARGET CGAL::NanoSVG_support)
  add_library(CGAL::NanoSVG_support INTERFACE IMPORTED)
  set_target_properties(CGAL::NanoSVG_support PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "CGAL_NANOSVG_ENABLED")
  target_link_libraries(CGAL::NanoSVG_support INTERFACE NanoSVG::nanosvg)
endif()
