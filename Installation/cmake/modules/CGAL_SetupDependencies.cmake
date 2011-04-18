include(CGAL_Macros)

message ( STATUS "External libraries supported: ${CGAL_SUPPORTING_3RD_PARTY_LIRARIES}")

foreach (lib ${CGAL_SUPPORTING_3RD_PARTY_LIRARIES})

  # 1) preconfiguration for non-lib CMakeLists.txt
 
  # TODO move "NOT LIB_FOUND" to second part of this loop 
  if (NOT ${lib}_FOUND OR WITH_${lib})
    preconfigure_lib( ${lib} )
  endif()


  # 2) handling for libraries required to build CGAL's libraries

  if ( ${lib} STREQUAL "GMPXX" )
  
    if ( MSVC AND NOT CGAL_AUTO_LINK_MPFR )
      add_to_cached_list(CGAL_3RD_PARTY_DEFINITIONS    -DCGAL_NO_AUTOLINK_MPFR   )
    endif()
    if ( MSVC AND NOT CGAL_AUTO_LINK_GMP )
      add_to_cached_list(CGAL_3RD_PARTY_DEFINITIONS    -DCGAL_NO_AUTOLINK_GMP   )
    endif()

    get_dependency_version(GMP)
    set( MPFR_DEPENDENCY_LIBRARIES   ${GMP_LIBRARIES} )
    set( MPFR_DEPENDENCY_INCLUDE_DIR ${GMP_INCLUDE_DIR} )
    get_dependency_version(MPFR)

    if ( NOT MPFR_FOUND )
      preconfigure_lib(MPFR)
    endif()

  endif()

  if ( ${lib} STREQUAL "LEDA" AND WITH_LEDA )
    include(CGAL_UseLEDA)
    uniquely_add_flags( CMAKE_CXX_FLAGS ${LEDA_CXX_FLAGS} )
  endif()

endforeach()

if( NOT GMP_FOUND )
  set(CGAL_NO_CORE ON)
  message( STATUS "CGAL_Core needs GMP, cannot be configured.")
endif( NOT GMP_FOUND )

message( STATUS "Preconfigured libraries: ${CGAL_3RD_PARTY_PRECONFIGURED}")

include(CGAL_SetupBoost)

if ( MSVC )
  add_to_cached_list(CGAL_3RD_PARTY_LIBRARIES "psapi.lib" )
endif( MSVC )



