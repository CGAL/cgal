if ( NOT CGAL_LEDA_SETUP )
  
  find_package( LEDA )
  
  if ( LEDA_FOUND )
  
    message( STATUS "LEDA include:     ${LEDA_INCLUDE_DIR}" )
    message( STATUS "LEDA libraries:   ${LEDA_LIBRARIES}" )
    message( STATUS "LEDA definitions: ${LEDA_DEFINITIONS}" )
    
    set ( CGAL_USE_LEDA 1 )
    
    include(CGAL_Macros)
    
    cache_set(CGAL_3RD_PARTY_INCLUDE_DIRS ${CGAL_3RD_PARTY_INCLUDE_DIRS} ${LEDA_INCLUDE_DIR} )
    cache_set(CGAL_3RD_PARTY_DEFINITIONS  ${CGAL_3RD_PARTY_DEFINITIONS}  ${LEDA_DEFINITIONS} )
    cache_set(CGAL_3RD_PARTY_LIBRARIES    ${CGAL_3RD_PARTY_LIBRARIES}    ${LEDA_LIBRARIES}   )
    
    uniquely_add_flags( CMAKE_CXX_FLAGS ${LEDA_CXX_FLAGS} )
    
    if ( BUILD_SHARED_LIBS )
      uniquely_add_flags( CMAKE_SHARED_LINKER_FLAGS ${LEDA_LINKER_FLAGS} )
    else()
      uniquely_add_flags( CMAKE_MODULE_LINKER_FLAGS ${LEDA_LINKER_FLAGS} )
    endif()
    
    uniquely_add_flags( CGAL_CXX_FLAGS ${LEDA_CXX_FLAGS} )
    
    get_dependency_version(LEDA)
    
  endif()
  
  
  set ( CGAL_LEDA_SETUP TRUE )
  
endif()

