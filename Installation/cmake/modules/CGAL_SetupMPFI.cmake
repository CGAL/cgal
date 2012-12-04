if ( NOT CGAL_MPFI_SETUP )
  
  find_package( MPFI )
  
  if ( MPFI_FOUND )
  
    message( STATUS "MPFI include:     ${MPFI_INCLUDE_DIR}" )
    message( STATUS "MPFI libraries:   ${MPFI_LIBRARIES}" )
    message( STATUS "MPFI definitions: ${MPFI_DEFINITIONS}" )
    
    set ( CGAL_USE_MPFI TRUE )
    
    include(CGAL_Macros)
    
    add_to_cached_list(CGAL_3RD_PARTY_INCLUDE_DIRS ${MPFI_INCLUDE_DIR} )
    add_to_cached_list(CGAL_3RD_PARTY_DEFINITIONS  ${MPFI_DEFINITIONS} )
    add_to_cached_list(CGAL_3RD_PARTY_LIBRARIES    ${MPFI_LIBRARIES}
                                                   ${MPFI_LINKER_FLAGS})

    uniquely_add_flags( CMAKE_CXX_FLAGS ${MPFI_CXX_FLAGS} )
    
  endif()
  
  
  set ( CGAL_MPFI_SETUP TRUE )
  
endif()
