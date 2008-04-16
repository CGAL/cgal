if ( NOT GENERATOR_SPECIFIC_SETTINGS_FILE_INCLUDED )
set( GENERATOR_SPECIFIC_SETTINGS_FILE_INCLUDED 1 )
 
  message( STATUS "Targetting ${CMAKE_GENERATOR}")
  
  if ( MSVC )
    # CGAL 3.2 level is 2
    SET(CMAKE_CXX_WARNING_LEVEL 2 CACHE STRING "C++ compiler warning level" FORCE)
    message( STATUS "Target build enviroment supports auto-linking" )
    set(AUTO_LINK_ENABLED TRUE)
  endif()

  if ( MSVC90 )  
    set(TOOLSET "vc90")
    message( STATUS "Using VC90 compiler." )
  elseif ( MSVC80 )  
    set(TOOLSET "vc80")
    message( STATUS "Using VC80 compiler." )
  elseif ( MSVC71 )
    set(TOOLSET "vc71")
    message( STATUS "Using VC71 compiler." )
  else()
    message( STATUS "Using ${CMAKE_CXX_COMPILER} compiler." )
  endif()

  if ( CMAKE_BUILD_TYPE )
    
    # If a build type is specified, VC project files contain only a configurartion for that build type
    set(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE} )
    message( STATUS "Build type: ${CMAKE_BUILD_TYPE}" )
    
  else()
  
    if ( CMAKE_CONFIGURATION_TYPES )
      message( STATUS "Build types: ${CMAKE_CONFIGURATION_TYPES}" )
    else()
      set(CMAKE_BUILD_TYPE Debug )
      message( STATUS "Build Type: ${CMAKE_BUILD_TYPE}" )
    endif()
  endif()
  
  if ( "${CMAKE_CONFIGURATION_TYPES}" MATCHES "Release" OR "${CMAKE_BUILD_TYPE}" MATCHES "Release" )
    set( CGAL_RELEASE TRUE )
  endif()
  
  if ( "${CMAKE_CONFIGURATION_TYPES}" MATCHES "Debug" OR "${CMAKE_BUILD_TYPE}" MATCHES "Debug" )
    set( CGAL_DEBUG TRUE )
  endif()
  
  if ( NOT "${CMAKE_CFG_INTDIR}" STREQUAL "." )
    set(HAS_CFG_INTDIR TRUE CACHE INTERNAL "Generator uses intermediate configuration directory" )
    message( STATUS "Generator uses intermediate configuration directory: ${CMAKE_CFG_INTDIR}" )
  endif()

  mark_as_advanced(CMAKE_CXX_WARNING_LEVEL)

endif()
