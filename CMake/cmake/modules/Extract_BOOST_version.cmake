include(ReadLines)
include(FindMatchingItem)

if ( Boost_FOUND )
  
  readlines("${Boost_INCLUDE_DIRS}/boost/version.hpp" BOOST_VERSION_FILE)
  
  if ( BOOST_VERSION_FILE )
  
    find_matching_item(BOOST_VERSION_FILE "BOOST_LIB_VERSION" BOOST_LIB_VERSION_LINE )
    
    string( REGEX MATCH "\".*\"$" BOOST_LIB_VERSION_STR2 ${BOOST_LIB_VERSION_LINE} )
    string( REPLACE "\"" "" BOOST_LIB_VERSION_STR1 ${BOOST_LIB_VERSION_STR2} )
    string( REPLACE "_" "." BOOST_VERSION ${BOOST_LIB_VERSION_STR1} )
    
  else()
  
    message( STATUS "WARNING: BOOST found but could not open ${Boost_INCLUDE_DIRS}/boost/version.hpp" )
    
    set ( BOOST_VERSION "unknown" )
    
  endif()
  
  message( STATUS "USING BOOST_VERSION = '${BOOST_VERSION}'" )
   
  
endif()

if ( BOOST_THREAD_FOUND )

  set ( BOOST_THREAD_VERSION ${BOOST_VERSION} )
  
  message( STATUS "USING BOOST_THREAD_VERSION = '${BOOST_THREAD_VERSION}'" )
  
endif()
