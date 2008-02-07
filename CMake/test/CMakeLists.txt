project(CGAL_TEST)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

file( GLOB list "*" )

foreach( entry ${list} )
  
  if ( IS_DIRECTORY ${entry} )
    
    if ( EXISTS ${entry}/CMakeLists.txt )
      message( STATUS "Configuring  ${entry} test" )
      add_subdirectory( ${entry} )
    else()
      message( STATUS "Skipping ${entry} test" )
    endif()
    
  endif()
  
endforeach()
