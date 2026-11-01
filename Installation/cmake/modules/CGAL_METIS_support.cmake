if(METIS_FOUND AND NOT TARGET CGAL::METIS_support)
  if(WIN32)
    find_package(GKlib)
    if (GKlib_FOUND AND (TARGET GKlib OR TARGET GKlib::GKlib))
      add_library(CGAL::METIS_support INTERFACE IMPORTED)
      set_target_properties(CGAL::METIS_support PROPERTIES
        INTERFACE_COMPILE_DEFINITIONS "CGAL_METIS_ENABLED"
        INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${METIS_LIBRARIES}")
      if (TARGET GKlib)
        target_link_libraries(CGAL::METIS_support INTERFACE GKlib)
      elseif(TARGET GKlib::GKlib)
        target_link_libraries(CGAL::METIS_support INTERFACE GKlib::GKlib)
      endif()
    else()
      message(STATUS "NOTICE: the METIS library requires the GKlib library where no compatible version has been found.")
    endif()
  else()
    add_library(CGAL::METIS_support INTERFACE IMPORTED)
    set_target_properties(CGAL::METIS_support PROPERTIES
      INTERFACE_COMPILE_DEFINITIONS "CGAL_METIS_ENABLED"
      INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "${METIS_LIBRARIES}")
  endif()
endif()
