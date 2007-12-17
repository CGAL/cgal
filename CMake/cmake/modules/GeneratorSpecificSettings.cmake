if ( NOT GENERATOR_SPECIFIC_SETTINGS_FILE_INCLUDED )
set( GENERATOR_SPECIFIC_SETTINGS_FILE_INCLUDED 1 )
 
  message( STATUS "Targetting ${CMAKE_GENERATOR}")
   
  if ( "${CMAKE_GENERATOR}" MATCHES "Visual Studio" )
    
    message( STATUS "Target build enviroment supports auto-linking" )
    set(AUTO_LINK_ENABLED TRUE)

  endif()

  # TOOLSET is used for mangling library names
  if ( "${CMAKE_GENERATOR}" MATCHES "2003" )
    set(TOOLSET "-vc71")
  elseif ( "${CMAKE_GENERATOR}" MATCHES "2005" )  
    set(TOOLSET "-vc80")
  elseif ( "${CMAKE_GENERATOR}" MATCHES "2008" )  
    set(TOOLSET "-vc90")
  endif()


  # Set Microsoft Visual C++ compilation warning level
  if (CMAKE_CXX_COMPILER MATCHES "^(CL|cl)")
    # CGAL 3.2 level is 2
    SET(CMAKE_CXX_WARNING_LEVEL 2	CACHE STRING "C++ compiler warning level" FORCE)
  endif()

  if(TOOLSET MATCHES "-vc" )
    set(COMMON_DEFINITIONS ${COMMON_DEFINITIONS} "-D_SECURE_SCL=0" )
    set(COMMON_DEFINITIONS ${COMMON_DEFINITIONS} "-D_CRT_SECURE_NO_DEPRECATE" )
    set(COMMON_DEFINITIONS ${COMMON_DEFINITIONS} "-DSCL_SECURE_NO_DEPRECATE" )
  endif()

  mark_as_advanced(CMAKE_CXX_WARNING_LEVEL)

endif()
