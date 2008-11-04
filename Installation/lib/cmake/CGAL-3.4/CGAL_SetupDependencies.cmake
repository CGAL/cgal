if ( MSVC )
  cache_set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} "psapi.lib" )
endif()

include(CGAL_SetupBoost)
include(CGAL_SetupGMP)
include(CGAL_SetupMPFR)
include(CGAL_SetupGMPXX)

option ( WITH_LEDA "Use the LEDA number types" ON )
if ( WITH_LEDA )
  include(CGAL_SetupLEDA)
endif()
