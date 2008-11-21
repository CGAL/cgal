if ( NOT CGAL_MPFR_SETUP )

  find_package( MPFR REQUIRED )
  
  message( STATUS "MPFR include:      ${MPFR_INCLUDE_DIR}" )
  message( STATUS "MPFR libraries:    ${MPFR_LIBRARIES}" )
  message( STATUS "MPFR definitions:  ${MPFR_DEFINITIONS}" )
  
  include(CGAL_Macros)
  
  set ( CGAL_USE_MPFR 1 )
  
  set( MPFR_DEPENDENCY_INCLUDE_DIR ${GMP_INCLUDE_DIR} )
  get_dependency_version(MPFR)
  
  add_to_list(CGAL_3RD_PARTY_INCLUDE_DIRS   PERSISTENT ${MPFR_INCLUDE_DIR}   )
  add_to_list(CGAL_3RD_PARTY_LIBRARIES_DIRS PERSISTENT ${MPFR_LIBRARIES_DIR} )
  add_to_list(CGAL_3RD_PARTY_DEFINITIONS    PERSISTENT ${MPFR_DEFINITIONS}   )
  
  if ( NOT MSVC )
    add_to_persistent_tagged_libraries(CGAL_3RD_PARTY_LIBRARIES ${MPFR_LIBRARIES} )
  endif()
  
  set( CGAL_MPFR_SETUP TRUE ) 
  
endif()


