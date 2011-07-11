project(CGAL_EXAMPLES)

cmake_minimum_required(VERSION 2.6.2)

foreach (package ${CGAL_CONFIGURED_PACKAGES})
  #message (STATUS "Current package: ${package}")
  file( GLOB listtmp "${package}/examples/*")
  list(APPEND list ${listtmp})
endforeach()

list( SORT list )

foreach( entry ${list} )
  
  if ( IS_DIRECTORY ${entry} )
    
    if ( EXISTS ${entry}/CMakeLists.txt )

      message( STATUS "Configuring examples in ${entry}" )

      string( REGEX REPLACE "${CMAKE_SOURCE_DIR}/*/examples/" "" DEMO_DIR "${entry}" )

      if( NOT "${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}") # out-of-source
       make_directory(${CMAKE_BINARY_DIR}/demo/${DEMO_DIR})
      endif()
 
      set(ADD_SUBDIR TRUE)
      if(EXISTS ${entry}/../../dont_submit)
        file(STRINGS ${entry}/../../dont_submit dont_submit_grep REGEX "^${DEMO_DIR}/?\$")
        if(dont_submit_grep) 
          set(ADD_SUBDIR FALSE)
        endif()
        file(STRINGS ${entry}/../../dont_submit dont_submit_grep REGEX "^examples/${DEMO_DIR}/?\$")
        if(dont_submit_grep) 
          set(ADD_SUBDIR FALSE)
        endif()
      endif()
      if(ADD_SUBDIR)
        add_subdirectory( ${entry} ${CMAKE_BINARY_DIR}/examples/${DEMO_DIR} )
      else()
        message(STATUS "examples/${DEMO_DIR} is in dont_submit")
      endif()

    endif()
    
  endif()
  
endforeach()

