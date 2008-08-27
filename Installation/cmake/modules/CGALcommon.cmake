# This allows else(), endif(), etc... (without repeating the expression)
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( "${CMAKE_SOURCE_DIR}" STREQUAL "${PROJECT_SOURCE_DIR}" )
  set( IS_TOP_LEVEL TRUE )
else()
  set( IS_TOP_LEVEL FALSE )
endif()  

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

  macro( cache_set var value )
    set ( ${var} ${value} CACHE INTERNAL "" FORCE )  
  endmacro()
  
  macro( cache_get var )
    set ( ${var} )  
  endmacro()
  
  macro( add_to_cached_list listname item )
    cache_get ( ${listname} )
    cache_set ( ${listname} "${${listname}};${item}" )
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
  
    if ( "${ARGC}" GREATER "1" )
      set( PKG ${ARGV1} )
    else()
      set( PKG ${LIB} )
    endif()
    
    if ( ${PKG}_FOUND )
    
      set ( ${LIB}_VERSION "unknown" )
      
      try_run( ${LIB}_RUN_RES
               ${LIB}_COMPILE_RES 
               ${CMAKE_BINARY_DIR} 
               ${CMAKE_SOURCE_DIR}/config/support/print_${LIB}_version.cpp 
               CMAKE_FLAGS -DINCLUDE_DIRECTORIES:STRING=${${PKG}_INCLUDE_DIR} 
                           -DLINK_LIBRARIES:STRING=${${PKG}_LIBRARIES}
                           -DLINK_DIRECTORIES:STRING=${${PKG}_LIBRARY_DIR}                           
               OUTPUT_VARIABLE ${LIB}_OUTPUT 
            )
            
      if ( ${LIB}_COMPILE_RES )
      
        if ( ${LIB}_RUN_RES EQUAL "0" )
        
          string( REGEX MATCH "version=.*\$" ${LIB}_VERSION_LINE ${${LIB}_OUTPUT}  )
          string( REPLACE "\n" "" ${LIB}_VERSION_LINE2 ${${LIB}_VERSION_LINE} )
          string( REPLACE "\r" "" ${LIB}_VERSION_LINE3 ${${LIB}_VERSION_LINE2} )
          string( REPLACE "version=" "" ${LIB}_VERSION ${${LIB}_VERSION_LINE3} )
          
        else()
        
          message( STATUS "WARNING: ${LIB} found but print_${LIB}_version.cpp exited with error condition: ${${LIB}_RUN_RES}" )
          message( STATUS "${PKG}_INCLUDE_DIR=${${PKG}_INCLUDE_DIR}" )
          message( STATUS "${PKG}_LIBRARIES=${${PKG}_LIBRARIES}" )
          message( STATUS "${PKG}_LIBRARY_DIR=${${PKG}_LIBRARY_DIR}" )
          message( STATUS "${${LIB}_OUTPUT}" )
          
        endif()
        
      else()
      
        message( STATUS "WARNING: ${LIB} found but could not compile print_${LIB}_version.cpp:")
        message( STATUS "${PKG}_INCLUDE_DIR=${${PKG}_INCLUDE_DIR}" )
        message( STATUS "${PKG}_LIBRARIES=${${PKG}_LIBRARIES}" )
        message( STATUS "${PKG}_LIBRARY_DIR=${${PKG}_LIBRARY_DIR}" )
        message( STATUS "${${LIB}_OUTPUT}" )
        
      endif() 
      
      message( STATUS "USING ${LIB}_VERSION = '${${LIB}_VERSION}'" )
  
    endif()
    
  endmacro()
  
  # CMAKE_ROOT must be properly configured, but is not by the CMake windows installer, so check here
  if (NOT CMAKE_ROOT)
    message( FATAL_ERROR "CMAKE_ROOT enviroment variable not set. It should point to the directory where CMake is installed.")
  endif()

  if ( COMMAND cmake_policy )
    cmake_policy( SET CMP0003 NEW )  
  endif()
  
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