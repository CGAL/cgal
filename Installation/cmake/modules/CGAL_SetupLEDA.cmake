if ( NOT LEDA_FOUND )
  
  find_package( LEDA )
  
  if ( LEDA_FOUND )
  
    message( STATUS "LEDA include:     ${LEDA_INCLUDE_DIR}" )
    message( STATUS "LEDA libraries:   ${LEDA_LIBRARIES}" )
    message( STATUS "LEDA definitions: ${LEDA_DEFINITIONS}" )
    
    set ( CGAL_USE_LEDA 1 )
    
    include(CGAL_Macros)
    
    include_directories( ${LEDA_INCLUDE_DIR} )
    
    link_libraries( ${LEDA_LIBRARIES} )
    
    add_definitions ( ${LEDA_DEFINITIONS} )
    
  endif()
  
endif()

