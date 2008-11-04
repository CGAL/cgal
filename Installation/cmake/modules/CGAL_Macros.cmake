if( NOT CGAL_MACROS_FILE_INCLUDED )
  set(CGAL_MACROS_FILE_INCLUDED 1 )
  
  macro(assert _arg )
    if ( NOT ${_arg} )
      message( FATAL_ERROR "Variable ${_arg} must be defined" ) 
    endif()
  endmacro()

  macro( hide_variable var )
    set ( ${var} ${${var}} CACHE INTERNAL "Variable hidden from user" FORCE )  
  endmacro()

  macro( cache_set var )
    if ( "${ARGC}" GREATER "1"  )
      set ( ${var} ${ARGN} CACHE INTERNAL "" FORCE )  
    endif()  
  endmacro()
  
  macro( cache_get var )
    set ( ${var} )  
  endmacro()
  
  macro( add_to_cached_list listname )
    cache_get ( ${listname} )
    cache_set ( ${listname} ${${listname}} ${ARGN} )
  endmacro()
  
  macro( at list idx var )
    list( LENGTH ${list} ${list}_length )
    if ( ${idx} LESS ${${list}_length} )
      list( GET ${list} ${idx} ${var} )
    else()
      set( ${var} "NOTFOUND" )    
    endif()
  endmacro()
  
  macro( found_in_list item_list item result )
    set( ${result} "FALSE" )
    foreach( element ${${item_list}} )
      if ( "${element}" STREQUAL "${item}" )
        set( ${result} "TRUE" )
      endif()
    endforeach()  
  endmacro()
  
  macro( uniquely_add_flags target_var )
    if ( "${ARGC}" GREATER "1"  )
      set( target_list "${${target_var}}" )
      separate_arguments( target_list )
      foreach( flag ${ARGN} )
        found_in_list( target_list ${flag} ${flag}_FOUND )
        if ( NOT ${flag}_FOUND )
          set( ${target_var} "${${target_var}} ${flag}" CACHE STRING "User-defined flags" FORCE )
        endif()  
      endforeach()
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
               CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${${PKG}_DEPENDENCY_INCLUDE_DIR};${${PKG}_INCLUDE_DIR}"
                           "-DLINK_LIBRARIES:STRING=${${PKG}_DEPENDENCY_LIBRARIES};${${PKG}_LIBRARIES}"
                           "-DLINK_DIRECTORIES:STRING=${${PKG}_DEPENDENCY_LIBRARY_DIR};${${PKG}_LIBRARY_DIR}"
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

  macro( create_CGALconfig_files )
    # FindCGAL and UseCGAL are platform specific so they are generated and stored in the binary folder.
    configure_file(${CGAL_MODULES_DIR}/CGALConfig_binary.cmake.in  ${CMAKE_BINARY_DIR}/CGALConfig.cmake       @ONLY IMMEDIATE)
    
    if ( SOURCE_INSTALL )
      configure_file(${CGAL_MODULES_DIR}/CGALConfig_install.cmake.source.in ${CMAKE_BINARY_DIR}/config/CGALConfig.cmake @ONLY IMMEDIATE)
    else()
      configure_file(${CGAL_MODULES_DIR}/CGALConfig_install.cmake.fhs.in    ${CMAKE_BINARY_DIR}/config/CGALConfig.cmake @ONLY IMMEDIATE)
    endif()
  endmacro()
  
  macro ( fetch_env_var VAR )
    if ( "${${VAR}}" STREQUAL "" )
      set( ${VAR}_env_value "$ENV{${VAR}}" )
      if ( NOT "${${VAR}_env_value}" STREQUAL "" )
        set( ${VAR} ${${VAR}_env_value} )
      endif()
    endif()
  endmacro()
  
endif()