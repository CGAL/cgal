project(CGAL_DEMOS)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

file( GLOB list "*" )

foreach( entry ${list} )
  
  if ( IS_DIRECTORY ${entry} )
    
    if ( EXISTS ${entry}/CMakeLists.txt )
      message( STATUS "Configuring  ${entry} demo" )
      add_subdirectory( ${entry} )
    endif()
    
  endif()
  
endforeach()
