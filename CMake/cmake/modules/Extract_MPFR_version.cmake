include(ReadLines)
include(FindMatchingItem)

if ( MPFR_FOUND )

  readlines(${MPFR_INCLUDE_DIR}/mpfr.h MPFR_H_FILE)
  
  find_matching_item(MPFR_H_FILE "MPFR_VERSION_STRING" MPFR_VERSION_STRING_LINE )
  
  string( REPLACE " " ";" MPFR_VERSION_STRING_LINE_LIST ${MPFR_VERSION_STRING_LINE} )
  
  list( GET MPFR_VERSION_STRING_LINE_LIST 2 MPFR_VERSION_STR )
  
  string( REPLACE "\"" "" MPFR_VERSION ${MPFR_VERSION_STR} )
  
  message( STATUS "USING MPFR_VERSION = '${MPFR_VERSION}'" )
  
endif()  

