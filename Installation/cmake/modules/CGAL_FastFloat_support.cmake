if(FastFloat_FOUND AND NOT TARGET CGAL::FastFloat_support)
  add_library(CGAL::FastFloat_support INTERFACE IMPORTED)
  set_target_properties(CGAL::FastFloat_support PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "CGAL_USE_FastFloat"
    INTERFACE_INCLUDE_DIRECTORIES "${FastFloat_INCLUDE_DIR}"
endif()
