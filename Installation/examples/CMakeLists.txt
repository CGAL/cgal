project(CGAL_EXAMPLES)

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

foreach (package ${CGAL_CONFIGURED_PACKAGES})
  #message (STATUS "Current package: ${package}")
  file( GLOB listtmp "${package}/examples/*")
  list(APPEND list ${listtmp})
endforeach()

list( SORT list )

foreach( entry ${list} )
  
  if ( IS_DIRECTORY ${entry} )
    
    if ( EXISTS ${entry}/CMakeLists.txt )
      message( STATUS "Configuring  ${entry} example" )
      add_subdirectory( ${entry} )
    endif()
    
  endif()
  
endforeach()

