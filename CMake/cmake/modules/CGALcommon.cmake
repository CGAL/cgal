# Common settings for CGAL cmake scripts
if( NOT CGAL_COMMON_FILE_INCLUDED )
  set(CGAL_COMMON_FILE_INCLUDED 1 )
  
  # This allows else(), endif(), etc... (without repeating the expression)
  set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

  macro(assert _arg )
    if ( NOT ${_arg} )
      message( FATAL_ERROR "Variable ${_arg} must be defined" ) 
    endif()
  endmacro()

  macro( hide_variable var )
    set ( ${var} ${${var}} CACHE INTERNAL "Variable hidden from user" FORCE )  
  endmacro()

  macro( at list idx var )
    list( LENGTH ${list} ${list}_length )
    if ( ${idx} LESS ${${list}_length} )
      list( GET ${list} ${idx} ${var} )
    else()
      set( ${var} "NOTFOUND" )    
    endif()
  endmacro()
  
  macro( get_dependency_version LIB )
  
    if ( ${ARG1} )
      set( PKG ${ARG1} )
    else()
      set( PKG ${LIB} )
    endif()
    
    if ( ${PKG}_FOUND )
    
      try_run( ${LIB}_RUN_RES
               ${LIB}_COMPILE_RES 
               ${CMAKE_BINARY_DIR} 
               ${CMAKE_SOURCE_DIR}/config/support/print_${LIB}_version.cpp 
               CMAKE_FLAGS -DINCLUDE_DIRECTORIES:STRING=${${PKG}_INCLUDE_DIR} 
                           -DLINK_LIBRARIES:STRING=${${PKG}_LIBRARIES} 
               OUTPUT_VARIABLE ${LIB}_OUTPUT 
            )
            
      if ( ${LIB}_COMPILE_RES )
      
        string( REGEX MATCH "version=.*\$" ${LIB}_VERSION_LINE ${${LIB}_OUTPUT}  )
        string( REPLACE "\n" "" ${LIB}_VERSION_LINE2 ${${LIB}_VERSION_LINE} )
        string( REPLACE "version=" "" ${LIB}_VERSION ${${LIB}_VERSION_LINE2} )
        
      else()
      
        message( STATUS "WARNING: ${LIB} found but could not execute print_${LIB}__version.cpp" )
    
        set ( ${LIB}_VERSION "unknown" )
        
      endif() 
      
      message( STATUS "USING ${LIB}_VERSION = '${${LIB}_VERSION}'" )
  
    endif()
    
  endmacro()
  
  # CMAKE_ROOT must be properly configured, but is not by the CMake windows installer, so check here
  if (NOT CMAKE_ROOT)
    message( FATAL_ERROR "CMAKE_ROOT enviroment variable not set. It should point to the directory where CMake is installed.")
  endif()

  CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5 FATAL_ERROR)

  if ( NOT BUILD_SHARED_LIBS )
    if ( WIN32 )
      set(BUILD_SHARED_LIBS OFF)
    else()
      set(BUILD_SHARED_LIBS ON)
    endif()
  endif()
  
  if ( BUILD_SHARED_LIBS )
    message( STATUS "Building shared libraries" )
  else()
    message( STATUS "Building static libraries" )
  endif()
  
  if ( WIN32 )
    find_program(CMAKE_UNAME uname /bin /usr/bin /usr/local/bin )
    if(CMAKE_UNAME)
      exec_program(uname ARGS -s OUTPUT_VARIABLE CMAKE_SYSTEM_NAME2)
      if ( CMAKE_SYSTEM_NAME2 MATCHES "CYGWIN" )
        message( STATUS "This is the Windows CMake running within the cygwin platform." )
        set( WIN32_CMAKE_ON_CYGWIN TRUE CACHE INTERNAL "This is the cygwin platform." )
      endif()
    endif()
    hide_variable(CMAKE_UNAME)
  endif()

  set(CMAKE_COLORMAKEFILE ON)
  
  # Needed by the testsuite results parser
  set(CMAKE_VERBOSE_MAKEFILE ON)

endif()