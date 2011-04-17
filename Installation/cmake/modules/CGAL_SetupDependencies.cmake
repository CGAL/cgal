include(CGAL_Macros)

# TODO sorted list of external libs

option( WITH_GMP "Use the GMP number types if available." ON )
if ( WITH_GMP )
  include(CGAL_SetupGMPXX)
  include(CGAL_SetupGMP)
endif( WITH_GMP )

if( NOT GMP_FOUND )
  set(CGAL_NO_CORE ON)
  message( STATUS "CGAL_Core needs GMP, cannot be configured.")
endif( NOT GMP_FOUND )

option ( WITH_LEDA "Use the LEDA number types if available." OFF )
if ( WITH_LEDA )
  include(CGAL_SetupLEDA)
endif( WITH_LEDA )

if ( WITH_MPFI )
  preconfigure_lib( MPFI )
endif( WITH_MPFI )

if ( WITH_RS )
  preconfigure_lib( RS )
endif( WITH_RS )

if ( WITH_NTL )
  preconfigure_lib( NTL )
endif( WITH_NTL )

message( STATUS "Preconfigured libraries: ${CGAL_3RD_PARTY_PRECONFIGURED}")

include(CGAL_SetupBoost)

if ( MSVC )
  add_to_cached_list(CGAL_3RD_PARTY_LIBRARIES "psapi.lib" )
endif( MSVC )



