option( WITH_GMP "Use the GMP number types if available." ON )

if ( WITH_GMP )
  include(CGAL_SetupGMPXX)
  include(CGAL_SetupGMP)
endif( WITH_GMP )

if( NOT GMP_FOUND )
  set(CGAL_NO_CORE ON)
endif( NOT GMP_FOUND )

option ( WITH_LEDA "Use the LEDA number types if available." OFF )
if ( WITH_LEDA )
  include(CGAL_SetupLEDA)
endif( WITH_LEDA )

include(CGAL_SetupBoost)

if ( MSVC )
  add_to_cached_list(CGAL_3RD_PARTY_LIBRARIES "psapi.lib" )
endif( MSVC )
