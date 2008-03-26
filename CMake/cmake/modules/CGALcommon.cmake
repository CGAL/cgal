# Common settings for CGAL cmake scripts
if( NOT CGAL_COMMON_FILE_INCLUDED )
  set(CGAL_COMMON_FILE_INCLUDED 1 )
  
  # This allows else(), endif(), etc... (without repeating the expression)
  set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

  # CMAKE_ROOT must be properly configured, but is not by the CMake windows installer, so check here
  if (NOT CMAKE_ROOT)
    message( FATAL_ERROR "CMAKE_ROOT enviroment variable not set. It should point to the directory where CMake is installed.")
  endif()

  # Check that the version of CMake is high enough.
  # CPack was introduced in cmake 2.4
  # FindQt3 is buggy in CMake 2.4.4.
  CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5 FATAL_ERROR)

  if ( NOT BUILD_SHARED_LIBS )
    if ( WIN32 )
      set(BUILD_SHARED_LIBS OFF)
    else()
      set(BUILD_SHARED_LIBS ON)
    endif()
  endif()
  
  if ( NOT CYGWIN )
    find_program(CMAKE_UNAME uname /bin /usr/bin /usr/local/bin )
    if(CMAKE_UNAME)
      exec_program(uname ARGS -s OUTPUT_VARIABLE CMAKE_SYSTEM_NAME2)
      if ( CMAKE_SYSTEM_NAME2 MATCHES "CYGWIN" )
        message( STATUS "This is the Windows CMake running within the cygwin platform." )
        set( CYGWIN_ON_WIN TRUE CACHE INTERNAL "This is the cygwin platform." )
      endif()
    endif()
  endif()

  # Just for fun
  set(CMAKE_COLORMAKEFILE ON)

  macro(assert _arg )
    if ( NOT ${_arg} )
      message( FATAL_ERROR "Variable ${_arg} must be defined" ) 
    endif()
  endmacro()

endif()