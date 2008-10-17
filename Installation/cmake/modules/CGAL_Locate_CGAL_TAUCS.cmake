if ( NOT CGAL_TAUCS_FOUND ) 
 
  if ( MSVC )
  
    if ( EXISTS "${CMAKE_SOURCE_DIR}/auxiliary/taucs" )
      set( CGAL_TAUCS_INCLUDE_DIR   "${CMAKE_SOURCE_DIR}/auxiliary/taucs/include")
      set( CGAL_TAUCS_LIBRARIES_DIR "${CMAKE_SOURCE_DIR}/auxiliary/taucs/lib"    )
      set( CGAL_TAUCS_FOUND TRUE )
      set_cache( TAUCS_IN_CGAL_AUXILIARY TRUE )
    endif()
    
  else()
  
    fetch_env_var(CGAL_TAUCS_DIR)
    
    if ( NOT "${CGAL_TAUCS_DIR}" STREQUAL "" )
      if ( EXISTS ${CGAL_TAUCS_DIR} )
        
        # Find out which platforms have been configured
        file( GLOB CGAL_TAUCS_BUILD_CONTENT_LIST RELATIVE "${CGAL_TAUCS_DIR}/build" "${CGAL_TAUCS_DIR}/build/*" )
        
        set( CGAL_TAUCS_PLATFORMS "" )
        
        foreach ( CGAL_TAUCS_BUILD_CONTENT ${CGAL_TAUCS_BUILD_CONTENT_LIST} )
          if ( IS_DIRECTORY "${CGAL_TAUCS_DIR}/build/${CGAL_TAUCS_BUILD_CONTENT}" )
            set( CGAL_TAUCS_PLATFORMS ${CGAL_TAUCS_PLATFORMS} ${CGAL_TAUCS_BUILD_CONTENT} )
          endif()
        endforeach()
 
        list( LENGTH CGAL_TAUCS_PLATFORMS CGAL_TAUCS_PLATFORMS_LEN )
        
        # For now we'll just pick up the first platform as the one to use.
        if ( CGAL_TAUCS_PLATFORMS_LEN GREATER "0" )
        
          list(GET CGAL_TAUCS_PLATFORMS 0 CGAL_TAUCS_PLATFORM )
          
          set( CGAL_TAUCS_INCLUDE_DIR   "${CGAL_TAUCS_DIR}/src" "${CGAL_TAUCS_DIR}/build/${CGAL_TAUCS_PLATFORM}")
          set( CGAL_TAUCS_LIBRARIES_DIR "${CGAL_TAUCS_DIR}/lib/${CGAL_TAUCS_PLATFORM}" )
          
          set( CGAL_TAUCS_FOUND TRUE )
          
        endif()
        
      endif()
      
    endif()  
    
  endif()
  
endif()

