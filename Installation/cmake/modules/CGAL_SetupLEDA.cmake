if ( NOT LEDA_FOUND )
  
  find_package( LEDA )
  
  if ( LEDA_FOUND )
  
    message( STATUS "LEDA include:     ${LEDA_INCLUDE_DIR}" )
    message( STATUS "LEDA libraries:   ${LEDA_LIBRARIES}" )
    message( STATUS "LEDA definitions: ${LEDA_DEFINITIONS}" )
    
    #get_dependency_version(LEDA)
    
    include_directories ( ${LEDA_INCLUDE_DIR} )     
          
    add_definitions( ${LEDA_DEFINITIONS} "-DCGAL_USE_LEDA=1" )
  
    link_directories( ${LEDA_LIBRARIES_DIR} )
    
    set( CGAL_3RD_PARTY_LIBRARIES  ${CGAL_3RD_PARTY_LIBRARIES} ${LEDA_LIBRARIES} )
    
  endif()
  
endif()

