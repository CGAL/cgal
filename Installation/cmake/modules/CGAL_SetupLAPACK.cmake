if ( NOT LAPACK_FOUND )
  
  find_package( LAPACK )
  
  if ( LAPACK_FOUND )
  
    message( STATUS "LAPACK libraries:   ${LAPACK_LIBRARIES}"   )
    message( STATUS "LAPACK definitions: ${LAPACK_DEFINITIONS}" )
    
    #get_dependency_version(LAPACK)
    
    add_definitions( ${LAPACK_DEFINITIONS} "-DCGAL_USE_LAPACK=1" )
  
    set( CGAL_3RD_PARTY_LIBRARIES  ${CGAL_3RD_PARTY_LIBRARIES} ${LAPACK_LIBRARIES} )
    
  endif()
  
endif()

