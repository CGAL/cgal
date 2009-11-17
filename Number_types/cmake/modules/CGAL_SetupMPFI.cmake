if( NOT CGAL_MPFI_SETUP )

  option( WITH_MPFI "Search for MPFI" ON )

  if( WITH_MPFI )

    find_package( MPFI )

    if( MPFI_FOUND )

      include(CGAL_Macros)

      message( STATUS "MPFI include:    ${MPFI_INCLUDE_DIR}" )
      message( STATUS "MPFI libraries:  ${MPFI_LIBRARIES}" )

      set ( CGAL_USE_MPFI 1 )
      add_definitions ( "-DCGAL_USE_MPFI" )

      add_to_cached_list(CGAL_3RD_PARTY_INCLUDE_DIRS   ${MPFI_INCLUDE_DIR}   )
      add_to_cached_list(CGAL_3RD_PARTY_LIBRARIES_DIRS ${MPFI_LIBRARIES_DIR} )

      if( NOT MSVC )
        add_to_cached_list(CGAL_3RD_PARTY_LIBRARIES ${MPFI_LIBRARIES} )
      endif( NOT MSVC )

    endif( MPFI_FOUND )

  endif( WITH_MPFI )

  set( CGAL_MPFI_SETUP TRUE )

endif( NOT CGAL_MPFI_SETUP )
