include(CGAL_Macros)

message ( STATUS "External libraries supported: ${CGAL_SUPPORTING_3RD_PARTY_LIRARIES}")

foreach (lib ${CGAL_SUPPORTING_3RD_PARTY_LIRARIES})

  #message (STATUS "Check ${lib}")

  # Part 1: Try to find lib

  # TODO-EBEB do all libs have a WITH_ by default?

  list(FIND CGAL_MANDATORY_3RD_PARTY_LIBRARIES "${lib}" POSITION)
  if ("${POSITION}" STRGREATER "-1" OR WITH_${lib})
   
    # message (STATUS "With ${lib} given")
    
    find_package( ${lib} )
   
    if (${lib}_FOUND) 
      message( STATUS "${lib} is preconfigured with use-file: ${${lib}_USE_FILE}") 
      message( STATUS "${lib} include:     ${${lib}_INCLUDE_DIR}" )
      message( STATUS "${lib} libraries:   ${${lib}_LIBRARIES}" )
      message( STATUS "${lib} definitions: ${${lib}_DEFINITIONS}" )
      message( STATUS "${lib} cxx flags:   ${${lib}_CXX_FLAGS}" )
   
      set ( CGAL_USE_${lib} 1 )

      # Part 2: Add some lib-specific definitions or obtain version
   
      if (${lib} STREQUAL "GMP") 
        if ( MSVC AND NOT CGAL_AUTO_LINK_GMP )
          add_to_cached_list(CGAL_3RD_PARTY_DEFINITIONS    -DCGAL_NO_AUTOLINK_GMP   )
        endif()
        get_dependency_version(GMP)
      endif()

      if (${lib} STREQUAL "MPFR") 
        if ( MSVC AND NOT CGAL_AUTO_LINK_MPFR )
          add_to_cached_list(CGAL_3RD_PARTY_DEFINITIONS    -DCGAL_NO_AUTOLINK_MPFR   )
        endif()
        set( MPFR_DEPENDENCY_INCLUDE_DIR ${GMP_INCLUDE_DIR} )
        set( MPFR_DEPENDENCY_LIBRARIES   ${GMP_LIBRARIES} )
        get_dependency_version(MPFR)
      endif()

      if (${lib} STREQUAL "LEDA") 
        # special case for LEDA - add a flag
        uniquely_add_flags( CMAKE_CXX_FLAGS ${LEDA_CXX_FLAGS} )
      endif()

    else() 
   
      if ("${POSITION}" STRGREATER "-1") 
        message( ERROR "CGAL requires ${lib} to be found" )
      endif()

    endif()

  endif()

endforeach()

if( NOT GMP_FOUND )
  set(CGAL_NO_CORE ON)
  message( STATUS "CGAL_Core needs GMP, cannot be configured.")
endif( NOT GMP_FOUND )

# finally setup Boost
include(CGAL_SetupBoost)

if ( MSVC )
  add_to_cached_list(CGAL_3RD_PARTY_LIBRARIES "psapi.lib" )
endif( MSVC )

