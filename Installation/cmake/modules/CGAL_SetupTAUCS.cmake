include(CGAL_SetupBLAS)
include(CGAL_SetupLAPACK)

if ( NOT TAUCS_FOUND )
  
  find_package( TAUCS )
  
  if ( TAUCS_FOUND )
  
    message( STATUS "TAUCS include:     ${TAUCS_INCLUDE_DIR}" )
    message( STATUS "TAUCS libraries:   ${TAUCS_LIBRARIES}" )
    message( STATUS "TAUCS definitions: ${TAUCS_DEFINITIONS}" )
    
    get_dependency_version(TAUCS)
    
    include_directories ( ${TAUCS_INCLUDE_DIR} )     
          
    add_definitions( ${TAUCS_DEFINITIONS} "-DCGAL_USE_TAUCS=1" )
  
    link_directories( ${TAUCS_LIBRARIES_DIR} )
    
    set( CGAL_3RD_PARTY_LIBRARIES  ${CGAL_3RD_PARTY_LIBRARIES} ${TAUCS_LIBRARIES} )
    
  endif()
  
endif()

