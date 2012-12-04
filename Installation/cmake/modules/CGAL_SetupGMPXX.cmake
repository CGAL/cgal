if ( NOT WIN32 AND NOT CGAL_GMPXX_SETUP )
  
  option( WITH_GMPXX "Search for GMPXX" ON )
  
  if ( WITH_GMPXX )
  
    find_package( GMPXX QUIET )
    
    if ( GMPXX_FOUND )
    
      message( STATUS "GMPXX include:    ${GMPXX_INCLUDE_DIR}" )
      message( STATUS "GMPXX libraries:  ${GMPXX_LIBRARIES}" )
      
      set ( CGAL_USE_GMPXX 1 )
      
      include(CGAL_Macros)
      
      add_to_cached_list(CGAL_3RD_PARTY_INCLUDE_DIRS ${GMPXX_INCLUDE_DIR} )
      add_to_cached_list(CGAL_3RD_PARTY_LIBRARIES    ${GMPXX_LIBRARIES}   )
      
    endif()
    
  endif()
  
  set( CGAL_GMPXX_SETUP TRUE )
  
endif()

