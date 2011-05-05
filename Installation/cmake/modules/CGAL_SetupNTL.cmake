if ( NOT CGAL_NTL_SETUP )
  
  find_package( NTL )
  
  if ( NTL_FOUND )
  
    message( STATUS "NTL include:     ${NTL_INCLUDE_DIR}" )
    message( STATUS "NTL libraries:   ${NTL_LIBRARIES}" )
    message( STATUS "NTL definitions: ${NTL_DEFINITIONS}" )
    
    set ( CGAL_USE_NTL TRUE )
    
    include(CGAL_Macros)
    
    add_to_cached_list(CGAL_3RD_PARTY_INCLUDE_DIRS ${NTL_INCLUDE_DIR} )
    add_to_cached_list(CGAL_3RD_PARTY_DEFINITIONS  ${NTL_DEFINITIONS} )
    add_to_cached_list(CGAL_3RD_PARTY_LIBRARIES    ${NTL_LIBRARIES}
                                                   ${NTL_LINKER_FLAGS})

    uniquely_add_flags( CMAKE_CXX_FLAGS ${NTL_CXX_FLAGS} )
    
  endif()
  
  
  set ( CGAL_NTL_SETUP TRUE )
  
endif()
