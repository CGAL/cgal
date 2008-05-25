include(ReadLines)
include(FindMatchingItem)
  
if ( ZLIB_FOUND )

  readlines(${ZLIB_INCLUDE_DIR}/zlib.h ZLIB_H_FILE)

  if ( ZLIB_H_FILE )
  
    find_matching_item(ZLIB_H_FILE "define ZLIB_VERSION" ZLIB_VERSION_LINE )
    
    string( REPLACE " " ";" ZLIB_VERSION_LINE_LIST ${ZLIB_VERSION_LINE} )

    at( ZLIB_VERSION_LINE_LIST 2 ZLIB_VERSION_STR )
    
    string( REPLACE "\"" "" ZLIB_VERSION ${ZLIB_VERSION_STR} )
    
  else()

    message( STATUS "WARNING: ZLIB found but could not open ${ZLIB_INCLUDE_DIR}/zlib.h" )
    
    set ( ZLIB_VERSION "unknown" )
    
  endif()
  
  message( STATUS "USING ZLIB_VERSION = '${ZLIB_VERSION}'" )
  
endif()

