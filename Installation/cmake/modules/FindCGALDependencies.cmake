if ( MSVC )
  cache_set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} "psapi.lib" )
endif()

include(CGAL_SetupBoost)
include(CGAL_SetupGMP)
include(CGAL_SetupMPFR)
include(CGAL_SetupGMPXX)
include(CGAL_SetupLEDA)
