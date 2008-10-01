macro( add_subdirectory_if cond dir)

  if ( ${cond} )
    message( STATUS "Configuring ${dir}. Set ${cond} to FALSE to unselect it." )
    add_subdirectory( ${dir} ${ARGN} )
  endif()

endmacro()

macro( optional_add_subdirectory dir option def)
  set( ${option}_ENV $ENV{${option}} )
  if ( NOT ${${option}_ENV} STREQUAL "" )
    message ( STATUS "${option}_ENV given as enviroment variable: ${${option}_ENV}" )
    set( ${option} ${${option}_ENV} CACHE BOOL "Select ${option}." FORCE )
  endif()
  option( ${option} "Select ${option}." ${def} )
  add_subdirectory_if( ${option} ${dir} ${ARGN} )
endmacro()

