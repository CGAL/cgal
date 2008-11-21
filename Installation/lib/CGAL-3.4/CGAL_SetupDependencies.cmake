include(CGAL_SetupGMPXX)
include(CGAL_SetupMPFR)
include(CGAL_SetupGMP)

option ( WITH_LEDA "Use the LEDA number types" ON )
if ( WITH_LEDA )
  include(CGAL_SetupLEDA)
endif()

include(CGAL_SetupBoost)

if ( MSVC )
  add_to_persistent_tagged_libraries(CGAL_3RD_PARTY_LIBRARIES "psapi.lib" )
endif()

