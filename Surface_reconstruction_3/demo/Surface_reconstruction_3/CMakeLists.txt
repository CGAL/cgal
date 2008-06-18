project(Surface_reconstruction_3)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

file( GLOB list "*" )

foreach( entry ${list} )

  if ( IS_DIRECTORY ${entry} )

    if ( EXISTS ${entry}/CMakeLists.txt )
      message( STATUS "Configuring  ${entry} demo" )
      add_subdirectory( ${entry} )
    else()
      message( STATUS "Skipping ${entry} demo" )
    endif()

  endif()

endforeach()
