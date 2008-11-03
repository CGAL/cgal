if ( NOT LEDA_FOUND )
  
  find_package( LEDA )
  
  if ( LEDA_FOUND )
  
    message( STATUS "LEDA include:     ${LEDA_INCLUDE_DIR}" )
    message( STATUS "LEDA libraries:   ${LEDA_LIBRARIES}" )
    message( STATUS "LEDA definitions: ${LEDA_DEFINITIONS}" )
    
    set ( CGAL_USE_LEDA 1 )
    
    include(CGAL_Macros)
    
    cache_set(CGAL_3RD_PARTY_INCLUDE_DIRS ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${LEDA_INCLUDE_DIR} )
    cache_set(CGAL_3RD_PARTY_DEFINITIONS  ${CGAL_3RD_PARTY_DEFINITIONS}  ${LEDA_DEFINITIONS} )
    cache_set(CGAL_3RD_PARTY_LIBRARIES    ${CGAL_3RD_PARTY_LIBRARIES}    ${LEDA_LIBRARIES}   )
    
    message( STATUS "USING LEDA_VERSION = '${CGAL_LEDA_VERSION}'" )
    
  endif()
  
endif()

