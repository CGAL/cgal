include(ReadLines)
include(FindMatchingItem)

if ( QT_FOUND )
  
  set( QT_VERSION ${qt_version_str} )
   
  message( STATUS "USING QT_VERSION = '${QT_VERSION}'" )
  
endif()

