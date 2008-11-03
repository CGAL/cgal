if ( NOT MPFR_FOUND )

  find_package( MPFR REQUIRED )
  
  message( STATUS "MPFR include:      ${MPFR_INCLUDE_DIR}" )
  message( STATUS "MPFR libraries:    ${MPFR_LIBRARIES}" )
  message( STATUS "MPFR definitions:  ${MPFR_DEFINITIONS}" )
  
  include(CGAL_Macros)
  
  set ( CGAL_USE_MPFR 1 )
  
  set( MPFR_DEPENDENCY_INCLUDE_DIR ${GMP_INCLUDE_DIR} )
  get_dependency_version(MPFR)
  
  cache_set(CGAL_3RD_PARTY_INCLUDE_DIRS   ${CGAL_3RD_PARTY_INCLUDE_DIRS}   ${MPFR_INCLUDE_DIR}   )
  cache_set(CGAL_3RD_PARTY_LIBRARIES_DIRS ${CGAL_3RD_PARTY_LIBRARIES_DIRS} ${MPFR_LIBRARIES_DIR} )
  cache_set(CGAL_3RD_PARTY_DEFINITIONS    ${CGAL_3RD_PARTY_DEFINITIONS}    ${MPFR_DEFINITIONS}   )
  
  if ( NOT MSVC )
    cache_set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} ${MPFR_LIBRARIES})
  endif()
  
endif()


