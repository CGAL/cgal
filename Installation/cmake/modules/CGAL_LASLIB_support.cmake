if(LASLIB_FOUND)
  if (NOT TARGET CGAL::LASLIB_support)
    if (NOT TARGET LASlib)
      # message(STATUS "Found using MODULE mode")
      add_library(CGAL::LASLIB_support INTERFACE IMPORTED)
      set_target_properties(CGAL::LASLIB_support PROPERTIES
        INTERFACE_COMPILE_DEFINITIONS "CGAL_LINKED_WITH_LASLIB"
        INTERFACE_INCLUDE_DIRECTORIES "${LASLIB_INCLUDE_DIR}"
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${LASLIB_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${LASLIB_LIBRARIES}")
    else()
      # message(STATUS "Found using CONFIG mode")
      add_library(CGAL::LASLIB_support INTERFACE IMPORTED)
      set_target_properties(CGAL::LASLIB_support  PROPERTIES
        INTERFACE_COMPILE_DEFINITIONS "CGAL_LINKED_WITH_LASLIB")
      target_link_libraries(CGAL::LASLIB_support INTERFACE LASlib)
    endif()
  endif()
endif()
