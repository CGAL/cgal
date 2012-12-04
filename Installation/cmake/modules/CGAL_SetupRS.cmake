if ( NOT CGAL_RS_SETUP )
  
  find_package( RS )
  
  if ( RS_FOUND )
  
    message( STATUS "RS include:     ${RS_INCLUDE_DIR}" )
    message( STATUS "RS libraries:   ${RS_LIBRARIES}" )
    message( STATUS "RS definitions: ${RS_DEFINITIONS}" )

    if( APPLE AND CMAKE_COMPILER_IS_GNUCXX )
      include( CGAL_VersionUtils )
      EXEC_PROGRAM( ${CMAKE_CXX_COMPILER}
                    ARGS -dumpversion
                    OUTPUT_VARIABLE RS_GXX_VERSION )
      VERSION_DECOMPOSE( ${RS_GXX_VERSION} GXX_MAJ GXX_MIN GXX_PAT GXX_TWE )
      IS_VERSION_LESS( "${GXX_MAJ}.${GXX_MIN}" "4.3" IS_OLD_GXX )
      if( IS_OLD_GXX )
        message( STATUS "TLS is not supported by g++<4.3 on Mac OS X" )
        add_to_cached_list(CGAL_3RD_PARTY_DEFINITIONS "-DCGAL_RS_NO_TLS" )
      endif( IS_OLD_GXX )
    endif( APPLE AND CMAKE_COMPILER_IS_GNUCXX )

    set ( CGAL_USE_RS 1 )
    
    include(CGAL_Macros)
    
    add_to_cached_list(CGAL_3RD_PARTY_INCLUDE_DIRS ${RS_INCLUDE_DIR} )
    add_to_cached_list(CGAL_3RD_PARTY_DEFINITIONS  ${RS_DEFINITIONS} "-DCGAL_USE_RS" )
    add_to_cached_list(CGAL_3RD_PARTY_LIBRARIES    ${RS_LIBRARIES}
                                                   ${RS_LINKER_FLAGS})

    uniquely_add_flags( CMAKE_CXX_FLAGS ${RS_CXX_FLAGS} )
    
  endif()

  # add rs3 parameters, if necessary (rs3 must be always after rsexport)
  if( RS3_FOUND )

    message( STATUS "RS3 include:       ${RS3_INCLUDE_DIR}" )
    message( STATUS "RS3 definitions:   ${RS3_DEFINITIONS}" )
    message( STATUS "RS3 libraries:     ${RS3_LIBRARIES}" )

    set ( CGAL_USE_RS3 1 )
    
    include(CGAL_Macros)
    
    add_to_cached_list(CGAL_3RD_PARTY_INCLUDE_DIRS ${RS3_INCLUDE_DIR} )
    add_to_cached_list(CGAL_3RD_PARTY_DEFINITIONS  ${RS3_DEFINITIONS} "-DCGAL_USE_RS3" )
    add_to_cached_list(CGAL_3RD_PARTY_LIBRARIES    ${RS3_LIBRARIES}
                                                   ${RS3_LINKER_FLAGS})

    uniquely_add_flags( CMAKE_CXX_FLAGS ${RS3_CXX_FLAGS} )
  
  endif( RS3_FOUND )


    
  
  
  set ( CGAL_RS_SETUP TRUE )
  
endif()
