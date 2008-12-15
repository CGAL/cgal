project(CGAL_DEMOS)

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

file( GLOB list "*" )

list( SORT list )

foreach( entry ${list} )
  
  if ( IS_DIRECTORY ${entry} )
    
    if ( EXISTS ${entry}/CMakeLists.txt )
      message( STATUS "Configuring  ${entry} demo" )
      add_subdirectory( ${entry} )
    endif()
    
  endif()
  
endforeach()
