include(MacroFindOptionalCGALDependency)

if ( MSVC )
  cache_set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} "psapi.lib" )
endif()

if ( NOT BUILD_SHARED_LIBS )
  set(Boost_USE_STATIC_LIBS ON)
endif()

set(Boost_FIND_VERSION 1.33.1 )
set(Boost_FIND_VERSION_MAJOR 1 )
set(Boost_FIND_VERSION_MINOR 33 )
set(Boost_FIND_VERSION_PATCH 1 )

find_package( Boost REQUIRED thread )

message( STATUS "Boost include:      ${Boost_INCLUDE_DIRS}" )
message( STATUS "Boost libraries:    ${Boost_LIBRARIES}" )
message( STATUS "Boost definitions:  ${Boost_DEFINITIONS}" )

cache_set(CGAL_3RD_PARTY_INCLUDE_DIRS  ${CGAL_3RD_PARTY_INCLUDE_DIRS}  ${Boost_INCLUDE_DIRS} "" )
cache_set(CGAL_3RD_PARTY_LIBRARIES_DIR ${CGAL_3RD_PARTY_LIBRARIES_DIR} ${Boost_LIBRARY_DIRS} "" )
cache_set(CGAL_3RD_PARTY_DEFINITIONS   ${CGAL_3RD_PARTY_DEFINITIONS}   ${Boost_DEFINITIONS}  "" )

if ( NOT AUTO_LINK_ENABLED )
  cache_set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES})
endif()

message( STATUS "USING BOOST_VERSION = '${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}'" )

find_optional_cgal_dependency(GMP  CGAL )
find_optional_cgal_dependency(MPFR CGAL )

if ( NOT WIN32 )
  find_optional_cgal_dependency(GMPXX CGAL )
endif()

if ( GMP_FOUND )
  get_dependency_version(GMP)
  message( STATUS "GMP include:      ${GMP_INCLUDE_DIR}" )
  message( STATUS "GMP libraries:    ${GMP_LIBRARIES}" )
  message( STATUS "GMP definitions:  ${GMP_DEFINITIONS}" )
endif()

if ( MPFR_FOUND )
  set( MPFR_DEPENDENCY_INCLUDE_DIR ${GMP_INCLUDE_DIR} )
  get_dependency_version(MPFR)
  message( STATUS "MPFR include:     ${MPFR_INCLUDE_DIR}" )
  message( STATUS "MPFR libraries:   ${MPFR_LIBRARIES}" )
  message( STATUS "MPFR definitions: ${MPFR_DEFINITIONS}" )
endif()