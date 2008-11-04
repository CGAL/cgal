if ( NOT Boost_FOUND )
  
  if ( NOT BUILD_SHARED_LIBS )
    set(Boost_USE_STATIC_LIBS ON)
  endif()
  
  set(Boost_FIND_VERSION 1.33.1 )
  set(Boost_FIND_VERSION_MAJOR 1 )
  set(Boost_FIND_VERSION_MINOR 33 )
  set(Boost_FIND_VERSION_PATCH 1 )
  
  find_package( Boost REQUIRED thread )
  
  message( STATUS "Boost include:     ${Boost_INCLUDE_DIRS}" )
  message( STATUS "Boost libraries:   ${Boost_LIBRARIES}" )
  message( STATUS "Boost definitions: ${Boost_DEFINITIONS}" )
  
  set ( CGAL_USE_BOOST 1 )
  
  include(CGAL_Macros)
  
  cache_set(CGAL_3RD_PARTY_INCLUDE_DIRS   ${CGAL_3RD_PARTY_INCLUDE_DIRS}   ${Boost_INCLUDE_DIRS} )
  cache_set(CGAL_3RD_PARTY_LIBRARIES_DIRS ${CGAL_3RD_PARTY_LIBRARIES_DIRS} ${Boost_LIBRARY_DIRS} )
  cache_set(CGAL_3RD_PARTY_DEFINITIONS    ${CGAL_3RD_PARTY_DEFINITIONS}    ${Boost_DEFINITIONS}  )
  
  if ( NOT MSVC )
    cache_set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES})
  endif()
  
  message( STATUS "USING BOOST_VERSION = '${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}'" )
  
endif()

