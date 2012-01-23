project(CGAL_EXAMPLES)

cmake_minimum_required(VERSION 2.6.2)

if (CGAL_BRANCH_BUILD) 

foreach (package ${CGAL_CONFIGURED_PACKAGES})
  #message (STATUS "Current package: ${package}")
  file( GLOB listtmp "${package}/examples/*")
  list(APPEND list ${listtmp})
endforeach()

else()

  file( GLOB list "*")

endif()

list( SORT list )

foreach( entry ${list} )

  if (NOT ${entry} MATCHES ".*\\.svn\$" AND IS_DIRECTORY ${entry} )

    process_CGAL_subdirectory("${entry}" examples example)
    # Note: process_CGAL_subdirectory is defined in ../CMakeLists.txt

  endif()
  
endforeach()
