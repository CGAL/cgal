add_subdirectory(CGAL)

create_CGALconfig_files()

set(CGAL_DIR ${CMAKE_BINARY_DIR})

macro( configure_component DIR COMPONENT )
  option( WITH_CGAL_${COMPONENT} "Enable CGAL component ${COMPONENT}" ON )
  if ( WITH_CGAL_${COMPONENT})
    add_subdirectory( ${DIR} )
  endif()
endmacro()

if(COMMAND unset)
  # Workaround: WITH_CGAL_Core can have been loaded from an old reference
  # cache Then remove it from the cache, and re-add it only if CGAL_Core
  # can be configured.
  unset(WITH_CGAL_Core CACHE)
endif()

if(NOT CGAL_NO_CORE)
  configure_component( CGALCore    Core    )
endif(NOT CGAL_NO_CORE)

configure_component( CGALimageIO ImageIO )
configure_component( CGALQt      Qt3     )
configure_component( CGALQt4     Qt4     )
