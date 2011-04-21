project(CGAL_DEMOS)

cmake_minimum_required(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

foreach (package ${CGAL_CONFIGURED_PACKAGES})
  #message (STATUS "Current package: ${package}")
  file( GLOB listtmp "${package}/demo/*")
  list(APPEND list ${listtmp})
endforeach()

list( SORT list )

foreach( entry ${list} )
  
  if ( IS_DIRECTORY ${entry} )
    
    if ( EXISTS ${entry}/CMakeLists.txt )

      message( STATUS "Configuring demo in ${entry}" )

      if( NOT "${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}") # out-of-source
        string( REGEX REPLACE "${CMAKE_SOURCE_DIR}/.*/demo/" "" DEMO_DIR "${entry}" )
        make_directory(${CMAKE_BINARY_DIR}/demo/${DEMO_DIR})
      endif()
 
      add_subdirectory( ${entry} ${CMAKE_BINARY_DIR}/demo/${DEMO_DIR} )

    Endif()
    
  endif()
  
endforeach()
