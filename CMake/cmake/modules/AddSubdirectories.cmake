macro ( add_subdirectories relpath )

  file( GLOB list "${relpath}/*" )
  
  foreach( entry ${list} )
    
    if ( IS_DIRECTORY ${entry} )
    
      add_directory( ${entry} )
      
    endif()
    
  endforeach()

endmacro()
