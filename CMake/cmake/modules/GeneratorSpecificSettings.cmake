if ( NOT GENERATOR_SPECIFIC_SETTINGS_FILE_INCLUDED )
set( GENERATOR_SPECIFIC_SETTINGS_FILE_INCLUDED 1 )
 
  message( STATUS "Targetting ${CMAKE_GENERATOR}")
   
  # TOOLSET is used for mangling library names
  
  if ( MSVC )
    # CGAL 3.2 level is 2
    SET(CMAKE_CXX_WARNING_LEVEL 2	CACHE STRING "C++ compiler warning level" FORCE)
    message( STATUS "Target build enviroment supports auto-linking" )
    set(AUTO_LINK_ENABLED TRUE)
  endif()
  
  if ( MSVC71 )
    set(TOOLSET "-vc71")
    message( STATUS "Using VC71 compiler." )
  elseif ( MSVC80 )  
    set(TOOLSET "-vc80")
    message( STATUS "Using VC80 compiler." )
  elseif ( MSVC90 )  
    set(TOOLSET "-vc90")
    message( STATUS "Using VC90 compiler." )
  else()
    message( STATUS "Using ${CMAKE_CXX_COMPILER} compiler." )
  endif()

  # Set Microsoft Visual C++ compilation warning level

  if( MSVC )
    
    if ( CMAKE_BUILD_TYPE )
      set(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE} )
    endif()
    
    message( STATUS "Build Types: ${CMAKE_CONFIGURATION_TYPES}" )
    
    set(COMMON_DEFINITIONS ${COMMON_DEFINITIONS} "-D_SECURE_SCL=0" )
    set(COMMON_DEFINITIONS ${COMMON_DEFINITIONS} "-D_CRT_SECURE_NO_DEPRECATE" )
    set(COMMON_DEFINITIONS ${COMMON_DEFINITIONS} "-DSCL_SECURE_NO_DEPRECATE" )
    
  else()  
  
    if(NOT CMAKE_BUILD_TYPE)
      set(CMAKE_BUILD_TYPE Release )
    endif()
    
    message( STATUS "Build Type: ${CMAKE_BUILD_TYPE}" )
    
  endif()

  if ( "${CMAKE_CFG_INTDIR}" MATCHES "/Debug|/Release|/RelWithDebInfo|/MinSizeRel" )
    set(HAS_CFG_INTDIR TRUE)
    message( STATUS "Generator uses intermediate configuration directory: ${CMAKE_CFG_INTDIR}" )
  endif()

  mark_as_advanced(CMAKE_CXX_WARNING_LEVEL)

endif()
