if(ortools_FOUND AND NOT TARGET CGAL::ortools_support)
  add_library(CGAL::ortools_support INTERFACE IMPORTED)
  set_target_properties(CGAL::ortools_support PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "CGAL_LINKED_WITH_ORTOOLS"
    INTERFACE_INCLUDE_DIRECTORIES "${ORTOOLS_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "ortools::ortools")
endif()
