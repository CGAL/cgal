# This module setups the compiler for the MPFI library.
# It assumes that find_package(MPFI) was already called.

if( NOT CGAL_MPFI_SETUP )

  if( MPFI_FOUND )

    message( STATUS "MPFI include:      ${MPFI_INCLUDE_DIR}" )
    message( STATUS "MPFI libraries:    ${MPFI_LIBRARIES}" )
    message( STATUS "MPFI definitions:  ${MPFI_DEFINITIONS}" )

    include_directories( ${MPFI_INCLUDE_DIR} )
    link_directories( ${MPFI_LIBRARIES_DIR} )
    add_definitions( ${MPFI_DEFINITIONS} "-DCGAL_USE_MPFI" )

    if( NOT MSVC )
      link_libraries( ${MPFI_LIBRARIES} )
    endif( NOT MSVC )

  endif( MPFI_FOUND )

  set( CGAL_MPFI_SETUP TRUE )

endif( NOT CGAL_MPFI_SETUP )
