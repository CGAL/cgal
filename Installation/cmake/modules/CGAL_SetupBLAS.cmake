if ( NOT BLAS_FOUND )
  
  find_package( BLAS )
  
  if ( BLAS_FOUND )
  
    message( STATUS "BLAS libraries:   ${BLAS_LIBRARIES}" )
    message( STATUS "BLAS definitions: ${BLAS_DEFINITIONS}" )
    
    #get_dependency_version(BLAS)
    
    add_definitions( ${BLAS_DEFINITIONS} "-DCGAL_USE_BLAS=1" )
  
    set( CGAL_3RD_PARTY_LIBRARIES  ${CGAL_3RD_PARTY_LIBRARIES} ${BLAS_LIBRARIES} )
    
  endif()
  
endif()

