macro ( add_subdirectories relpath type )

  file( GLOB list "${relpath}/*" )
  
  foreach( entry ${list} )
    
    if ( IS_DIRECTORY ${entry} )
      
      if ( EXISTS ${entry}/CMakeLists.txt )
        message( STATUS "Configuring  ${entry} ${type}" )
        add_subdirectory( ${entry} )
      else()
        message( STATUS "Skipping ${entry} ${type}" )
      endif()
      
    endif()
    
  endforeach()

endmacro()
