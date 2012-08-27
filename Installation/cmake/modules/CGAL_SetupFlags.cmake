if ( NOT CGAL_SETUP_FLAGS_INCLUDED )
  set( CGAL_SETUP_FLAGS_INCLUDED 1 )

#
# Set CGAL_DONT_OVERRIDE_CMAKE_FLAGS to TRUE if you need to prevent the CGAL configuration to
# override the flags used to build the libraries
#
set( CGAL_DONT_OVERRIDE_CMAKE_FLAGS_DESCRIPTION
     "Set this to TRUE if you want to define or modify any of CMAKE_*_FLAGS. When this is FALSE, all the CMAKE_*_FLAGS flags are overriden with the values used when building the CGAL libs. For CGAL_*_flags (used for ADDITIONAL flags) , there is no need to set this to TRUE." 
   )

option( CGAL_DONT_OVERRIDE_CMAKE_FLAGS 
        ${CGAL_DONT_OVERRIDE_CMAKE_FLAGS_DESCRIPTION}
        FALSE 
      )

if ( CGAL_CONFIG_LOADED AND NOT CGAL_DONT_OVERRIDE_CMAKE_FLAGS )

  typed_cache_set ( STRING "Build type: Release or Debug" CMAKE_BUILD_TYPE "${CGAL_BUILD_TYPE_INIT}" )

  string( TOUPPER "${CMAKE_BUILD_TYPE}" CGAL_BUILD_TYPE_UPPER )
  
  if ( CGAL_BUILD_SHARED_LIBS )
    set( CGAL_LINKER_FLAGS_TYPE SHARED )
  else()
    set( CGAL_LINKER_FLAGS_TYPE MODULE )
  endif()
  
  typed_cache_set ( STRING "C++ compiler flags for both Release and Debug"   CMAKE_CXX_FLAGS                                 "${CGAL_CXX_FLAGS_INIT}"                                                       )
  typed_cache_set ( STRING "C++ compiler flags for ${CGAL_BUILD_TYPE_UPPER}" CMAKE_CXX_FLAGS_${CGAL_BUILD_TYPE_UPPER}        "${CGAL_CXX_FLAGS_${CGAL_BUILD_TYPE_UPPER}_INIT}"                              )
  typed_cache_set ( STRING "Linker flags for both Release and Debug"         CMAKE_EXE_LINKER_FLAGS                          "${CGAL_${CGAL_LINKER_FLAGS_TYPE}_LINKER_FLAGS_INIT}"                          )
  typed_cache_set ( STRING "Linker flags for ${CGAL_BUILD_TYPE_UPPER}"       CMAKE_EXE_LINKER_FLAGS_${CGAL_BUILD_TYPE_UPPER} "${CGAL_${CGAL_LINKER_FLAGS_TYPE}_LINKER_FLAGS_${CGAL_BUILD_TYPE_UPPER}_INIT}" )
  
endif()

typed_cache_set( BOOL ${CGAL_DONT_OVERRIDE_CMAKE_FLAGS_DESCRIPTION} CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE )

uniquely_add_flags( CMAKE_CXX_FLAGS                   ${CGAL_CXX_FLAGS}                   )
uniquely_add_flags( CMAKE_CXX_FLAGS_RELEASE           ${CGAL_CXX_FLAGS_RELEASE}           )
uniquely_add_flags( CMAKE_CXX_FLAGS_DEBUG             ${CGAL_CXX_FLAGS_DEBUG}             )
uniquely_add_flags( CMAKE_MODULE_LINKER_FLAGS         ${CGAL_MODULE_LINKER_FLAGS}         )
uniquely_add_flags( CMAKE_MODULE_LINKER_FLAGS_RELEASE ${CGAL_MODULE_LINKER_FLAGS_RELEASE} )
uniquely_add_flags( CMAKE_MODULE_LINKER_FLAGS_DEBUG   ${CGAL_MODULE_LINKER_FLAGS_DEBUG}   )
uniquely_add_flags( CMAKE_SHARED_LINKER_FLAGS         ${CGAL_SHARED_LINKER_FLAGS}         )
uniquely_add_flags( CMAKE_SHARED_LINKER_FLAGS_RELEASE ${CGAL_SHARED_LINKER_FLAGS_RELEASE} )
uniquely_add_flags( CMAKE_SHARED_LINKER_FLAGS_DEBUG   ${CGAL_SHARED_LINKER_FLAGS_DEBUG}   )
uniquely_add_flags( CMAKE_EXE_LINKER_FLAGS            ${CGAL_EXE_LINKER_FLAGS}            )
uniquely_add_flags( CMAKE_EXE_LINKER_FLAGS_RELEASE    ${CGAL_EXE_LINKER_FLAGS_RELEASE}    )
uniquely_add_flags( CMAKE_EXE_LINKER_FLAGS_DEBUG      ${CGAL_EXE_LINKER_FLAGS_DEBUG}      )

# Set a default build type if none is given
if ( NOT CMAKE_BUILD_TYPE )
  if( RUNNING_CGAL_AUTO_TEST )
    typed_cache_set ( STRING "Build type: Release or Debug" CMAKE_BUILD_TYPE Debug   )
  else ()
    typed_cache_set ( STRING "Build type: Release or Debug" CMAKE_BUILD_TYPE Release )
  endif()
endif()

if( RUNNING_CGAL_AUTO_TEST )
  add_definitions(-DCGAL_TEST_SUITE)
endif()

if ( NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Release" AND NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )
  message( FATAL_ERROR "${CMAKE_BUILD_TYPE} is not a valid build type: only Release or Debug is allowed" )
endif()

message( STATUS "Build type: ${CMAKE_BUILD_TYPE}" )

string( TOUPPER "${CMAKE_BUILD_TYPE}" CGAL_BUILD_TYPE_UPPER )

message( STATUS "USING CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CGAL_BUILD_TYPE_UPPER}}'" )

if ( CGAL_BUILDING_LIBS )
  if ( BUILD_SHARED_LIBS )
    message( STATUS "USING LDFLAGS = '${CMAKE_SHARED_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS_${CGAL_BUILD_TYPE_UPPER}}'" )
  else()
    message( STATUS "USING LDFLAGS = '${CMAKE_STATIC_LINKER_FLAGS} ${CMAKE_STATIC_LINKER_FLAGS_${CGAL_BUILD_TYPE_UPPER}}'" )
  endif()
else()
  message( STATUS "USING EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_${CGAL_BUILD_TYPE_UPPER}}'" )
endif()

endif()

