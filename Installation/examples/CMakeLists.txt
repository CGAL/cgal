project(CGAL_EXAMPLES)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

file( GLOB list "*" )

list( SORT list )

foreach( entry ${list} )
  
  if ( IS_DIRECTORY ${entry} )
    
    if ( EXISTS ${entry}/CMakeLists.txt )
      message( STATUS "Configuring  ${entry} example" )
      add_subdirectory( ${entry} )
    endif()
    
  endif()
  
endforeach()

