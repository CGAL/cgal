if ( NOT WIN32 AND NOT GMPXX_FOUND )
  
  option( WITH_GMPXX "Search for GMPXX" ON )
  
  if ( WITH_GMPXX )
  
    find_package( GMPXX QUIET )
    
    if ( GMPXX_FOUND )
    
      message( STATUS "GMPXX include:    ${GMPXX_INCLUDE_DIR}" )
      message( STATUS "GMPXX libraries:  ${GMPXX_LIBRARIES}" )
      
      set ( CGAL_USE_GMPXX 1 )
      
      include(CGAL_Macros)
      
      cache_set(CGAL_3RD_PARTY_INCLUDE_DIRS   ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${GMPXX_INCLUDE_DIR}  "" )
      
      if ( NOT MSVC )
        cache_set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} ${GMPXX_LIBRARIES})
      endif()
      
    endif()
    
  endif()
  
endif()

