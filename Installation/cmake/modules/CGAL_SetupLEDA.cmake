if ( NOT CGAL_LEDA_SETUP )
  
  find_package( LEDA )
  
  if ( LEDA_FOUND )
  
    message( STATUS "LEDA include:     ${LEDA_INCLUDE_DIR}" )
    message( STATUS "LEDA libraries:   ${LEDA_LIBRARIES}" )
    message( STATUS "LEDA definitions: ${LEDA_DEFINITIONS}" )
    
    set ( CGAL_USE_LEDA 1 )
    
    include(CGAL_Macros)
    
    add_to_cached_list(CGAL_3RD_PARTY_INCLUDE_DIRS ${LEDA_INCLUDE_DIR} )
    add_to_cached_list(CGAL_3RD_PARTY_DEFINITIONS  ${LEDA_DEFINITIONS} )
    add_to_cached_list(CGAL_3RD_PARTY_LIBRARIES    ${LEDA_LIBRARIES}
                                                   ${LEDA_LINKER_FLAGS})

    uniquely_add_flags( CMAKE_CXX_FLAGS ${LEDA_CXX_FLAGS} )
    
  endif()
  
  
  set ( CGAL_LEDA_SETUP TRUE )
  
endif()
