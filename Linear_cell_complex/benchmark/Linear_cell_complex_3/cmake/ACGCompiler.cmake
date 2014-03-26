################################################################################
# Custom settings for compiler flags and similar
################################################################################

if (UNIX)

  set ( ADDITIONAL_CXX_DEBUG_FLAGS )
  set ( ADDITIONAL_CXX_RELEASE_FLAGS )
  set ( ADDITIONAL_CXX_RELWITHDEBINFO_FLAGS )

  set ( ADDITIONAL_C_DEBUG_FLAGS )
  set ( ADDITIONAL_C_RELEASE_FLAGS )
  set ( ADDITIONAL_C_RELWITHDEBINFO_FLAGS )

  ################################################################################
  # Defaults
  ################################################################################

  # add our standard flags for Template inclusion
  list(APPEND ADDITIONAL_CXX_DEBUG_FLAGS          "-DINCLUDE_TEMPLATES" )
  list(APPEND ADDITIONAL_CXX_RELEASE_FLAGS        "-DINCLUDE_TEMPLATES" )
  list(APPEND ADDITIONAL_CXX_RELWITHDEBINFO_FLAGS "-DINCLUDE_TEMPLATES" )

  # add our standard flags for Template inclusion
  list(APPEND ADDITIONAL_C_DEBUG_FLAGS            "-DINCLUDE_TEMPLATES" )
  list(APPEND ADDITIONAL_C_RELEASE_FLAGS          "-DINCLUDE_TEMPLATES" )
  list(APPEND ADDITIONAL_C_RELWITHDEBINFO_FLAGS   "-DINCLUDE_TEMPLATES" )

  # Increase the template depth as this might be exceeded from time to time
  IF( NOT CMAKE_SYSTEM MATCHES "SunOS*")
    list(APPEND ADDITIONAL_CXX_DEBUG_FLAGS          "-ftemplate-depth-100" )
    list(APPEND ADDITIONAL_CXX_RELEASE_FLAGS        "-ftemplate-depth-100" )
    list(APPEND ADDITIONAL_CXX_RELWITHDEBINFO_FLAGS "-ftemplate-depth-100" )
  ENDIF()

  ################################################################################
  # OS Defines
  ################################################################################

  if (APPLE)
    add_definitions( -DARCH_DARWIN )
  endif()

  ################################################################################
  # Build/Release Defines
  ################################################################################
  IF( NOT CMAKE_SYSTEM MATCHES "SunOS*")
    list(APPEND ADDITIONAL_CXX_DEBUG_FLAGS          "-DDEBUG" )
    list(APPEND ADDITIONAL_CXX_RELEASE_FLAGS        "-DNDEBUG" )
    list(APPEND ADDITIONAL_CXX_RELWITHDEBINFO_FLAGS "-DDEBUG" )

    list(APPEND ADDITIONAL_C_DEBUG_FLAGS            "-DDEBUG" )
    list(APPEND ADDITIONAL_C_RELEASE_FLAGS          "-DNDEBUG" )
    list(APPEND ADDITIONAL_C_RELWITHDEBINFO_FLAGS   "-DDEBUG" )
  ENDIF()

  ################################################################################
  # Warnings
  ################################################################################

  IF( NOT CMAKE_SYSTEM MATCHES "SunOS*")
    list(APPEND ADDITIONAL_CXX_DEBUG_FLAGS          "-W" "-Wall" "-Wno-unused" )
    list(APPEND ADDITIONAL_CXX_RELEASE_FLAGS        "-W" "-Wall" "-Wno-unused" )
    list(APPEND ADDITIONAL_CXX_RELWITHDEBINFO_FLAGS "-W" "-Wall" "-Wno-unused" )

    list(APPEND ADDITIONAL_C_DEBUG_FLAGS            "-W" "-Wall" "-Wno-unused" )
    list(APPEND ADDITIONAL_C_RELEASE_FLAGS          "-W" "-Wall" "-Wno-unused" )
    list(APPEND ADDITIONAL_C_RELWITHDEBINFO_FLAGS   "-W" "-Wall" "-Wno-unused" )
  ENDIF()

  if (APPLE)
    list(APPEND ADDITIONAL_CXX_DEBUG_FLAGS          "-Wno-non-virtual-dtor" )
    list(APPEND ADDITIONAL_CXX_RELEASE_FLAGS        "-Wno-non-virtual-dtor" )
    list(APPEND ADDITIONAL_CXX_RELWITHDEBINFO_FLAGS "-Wno-non-virtual-dtor" )
  endif ()

  ################################################################################
  # STL Vector checks
  ################################################################################

  # Pre initialize stl vector check variable
  if ( NOT STL_VECTOR_CHECKS )
    set ( STL_VECTOR_CHECKS false CACHE BOOL "Include full stl vector checks in debug mode (This option is only used in debug Mode!)" )
  endif ( NOT STL_VECTOR_CHECKS )

  # Add a flag to check stl vectors in debugging mode
  if ( STL_VECTOR_CHECKS AND NOT CMAKE_SYSTEM MATCHES "SunOS*"  )
    list(APPEND ADDITIONAL_CXX_DEBUG_FLAGS          "-D_GLIBCXX_DEBUG" )
    list(APPEND ADDITIONAL_CXX_DEBUG_FLAGS          "-D_GLIBCXX_DEBUG_PEDANTIC")
    list(APPEND ADDITIONAL_CXX_RELWITHDEBINFO_FLAGS "-D_GLIBCXX_DEBUG" )
    list(APPEND ADDITIONAL_CXX_RELWITHDEBINFO_FLAGS "-D_GLIBCXX_DEBUG_PEDANTIC")

    list(APPEND ADDITIONAL_C_DEBUG_FLAGS            "-D_GLIBCXX_DEBUG" )
    list(APPEND ADDITIONAL_C_DEBUG_FLAGS            "-D_GLIBCXX_DEBUG_PEDANTIC")
    list(APPEND ADDITIONAL_C_RELWITHDEBINFO_FLAGS   "-D_GLIBCXX_DEBUG" )
    list(APPEND ADDITIONAL_C_RELWITHDEBINFO_FLAGS   "-D_GLIBCXX_DEBUG_PEDANTIC")
  endif()

  ################################################################################
  # Process the additional flags:
  ################################################################################

  # Add the debug flags
  foreach( flag ${ADDITIONAL_CXX_DEBUG_FLAGS} )
    if( NOT CMAKE_CXX_FLAGS_DEBUG MATCHES "${flag}" )
      set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${flag} ")
    endif()
  endforeach()

  # Add the release flags
  foreach( flag ${ADDITIONAL_CXX_RELEASE_FLAGS} )
    if( NOT CMAKE_CXX_FLAGS_RELEASE MATCHES "${flag}" )
      set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${flag} ")
    endif()
  endforeach()

  # Add the release with debug info flags
  foreach( flag ${ADDITIONAL_CXX_RELWITHDEBINFO_FLAGS} )
    if( NOT CMAKE_CXX_FLAGS_RELWITHDEBINFO MATCHES "${flag}" )
      set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${flag} ")
    endif()
  endforeach()

  # Add the debug flags
  foreach( flag ${ADDITIONAL_C_DEBUG_FLAGS} )
    if( NOT CMAKE_C_FLAGS_DEBUG MATCHES "${flag}" )
      set( CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${flag} ")
    endif()
  endforeach()

  # Add the release flags
  foreach( flag ${ADDITIONAL_C_RELEASE_FLAGS} )
    if( NOT CMAKE_C_FLAGS_RELEASE MATCHES "${flag}" )
      set( CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${flag} ")
    endif()
  endforeach()

  # Add the release with debug info flags
  foreach( flag ${ADDITIONAL_C_RELWITHDEBINFO_FLAGS} )
    if( NOT CMAKE_C_FLAGS_RELWITHDEBINFO MATCHES "${flag}" )
      set( CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${flag} ")
    endif()
  endforeach()

  #TODO : Test and remove it?!
  IF( CMAKE_SYSTEM MATCHES "SunOS*")
    set (CMAKE_CFLAGS_RELEASE "-xO3")
    set (CMAKE_CXX_FLAGS_RELEASE "-xO3")
  endif ( CMAKE_SYSTEM MATCHES "SunOS*" )

endif ()
