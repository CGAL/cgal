macro( add_subdirectory_if cond dir)

  if ( ${cond} )
    message( STATUS "Configuring ${dir}. Set ${cond} to FALSE to unselect it." )
    add_subdirectory( ${dir} ${ARGN} )
  else()  
    message( STATUS "Skipping ${dir}. Set ${cond} to TRUE to select it." )
  endif()

endmacro()

macro( optional_add_subdirectory dir def)
  option( WITH_${dir} "Select ${dir} package." ${def} )
  add_subdirectory_if( WITH_${dir} ${dir} ${ARGN} )
endmacro()

