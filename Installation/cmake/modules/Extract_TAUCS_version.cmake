include(ReadLines)
include(FindMatchingItem)
  
if ( TAUCS_FOUND )
  
  # TAUCS provides no version number :-(
  #/ Version 1 is obsolete, thus we assume version 2 (latest is 2.2 on 03/2006)
  set( TAUCS_VERSION "2.x" )
  
  message( STATUS "USING TAUCS_VERSION = '${TAUCS_VERSION}'" )
  
endif()

