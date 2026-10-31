if(NOT TARGET CGAL::Scotch_support)
  # First try config mode in case Scotch was built and installed with CMake
  find_package(Scotch CONFIG QUIET)

  if(TARGET Scotch::scotchmetis OR TARGET Scotch::scotch)
    add_library(CGAL::Scotch_support INTERFACE IMPORTED)
    set_target_properties(CGAL::Scotch_support PROPERTIES
      INTERFACE_COMPILE_DEFINITIONS "CGAL_SCOTCH_ENABLED")
    if(TARGET Scotch::scotchmetis)
      target_link_libraries(CGAL::Scotch_support INTERFACE Scotch::scotchmetis)
    else()
      target_link_libraries(CGAL::Scotch_support INTERFACE Scotch::scotch)
    endif()
  else()
    # Otherwise look for system installation paths (e.g. Debian/Ubuntu packages)
    find_path(SCOTCH_INCLUDE_DIR NAMES scotch.h PATH_SUFFIXES scotch)
    find_path(SCOTCHMETIS_INCLUDE_DIR NAMES metis.h PATH_SUFFIXES metis)
    find_library(SCOTCH_LIBRARY NAMES scotch)
    find_library(SCOTCHERR_LIBRARY NAMES scotcherr)
    find_library(SCOTCHMETIS_LIBRARY NAMES scotchmetisv5 scotchmetis)

    if(SCOTCH_INCLUDE_DIR AND SCOTCHMETIS_INCLUDE_DIR AND SCOTCH_LIBRARY AND SCOTCHMETIS_LIBRARY)
      set(Scotch_FOUND TRUE CACHE INTERNAL "Scotch library found")
      set(SCOTCH_LIBRARIES ${SCOTCHMETIS_LIBRARY} ${SCOTCH_LIBRARY} ${SCOTCHERR_LIBRARY} CACHE INTERNAL "Scotch libraries")
      set(SCOTCH_INCLUDE_DIRS ${SCOTCHMETIS_INCLUDE_DIR} ${SCOTCH_INCLUDE_DIR} CACHE INTERNAL "Scotch include directories")

      add_library(CGAL::Scotch_support INTERFACE IMPORTED)
      set_target_properties(CGAL::Scotch_support PROPERTIES
        INTERFACE_COMPILE_DEFINITIONS "CGAL_SCOTCH_ENABLED"
        INTERFACE_COMPILE_OPTIONS "-include;scotch.h"
        INTERFACE_INCLUDE_DIRECTORIES "${SCOTCH_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${SCOTCH_LIBRARIES}")
    endif()
  endif()
endif()
