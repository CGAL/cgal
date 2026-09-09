if ( NOT CGAL_GENERATOR_SPECIFIC_SETTINGS_FILE_INCLUDED )
  set( CGAL_GENERATOR_SPECIFIC_SETTINGS_FILE_INCLUDED 1 )

  message( STATUS "Targeting ${CMAKE_GENERATOR}")

  if ( MSVC )
    message( STATUS "Target build environment supports auto-linking" )
    set(CGAL_AUTO_LINK_ENABLED TRUE)
  endif()

  if ( MSVC_TOOLSET_VERSION )
    set(CGAL_TOOLSET "vc${MSVC_TOOLSET_VERSION}")
    message( STATUS "Using VC toolset ${MSVC_TOOLSET_VERSION}." )
  elseif ( MSVC )
    message( AUTHOR_WARNING
      "MSVC detected but MSVC_TOOLSET_VERSION is not set. "
      "CGAL_TOOLSET will not be configured. "
      "This may indicate an unsupported CMake version or toolchain." )
  else()
    message( STATUS "Using ${CMAKE_CXX_COMPILER} compiler." )
  endif()


  # From james Bigler, in the cmake users list.
  IF (APPLE)
    execute_process(COMMAND uname -v
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    OUTPUT_VARIABLE DARWIN_VERSION)
    string(REGEX MATCH "[0-9]+" DARWIN_VERSION ${DARWIN_VERSION})
    message(STATUS "Running in macOS DARWIN_VERSION=${DARWIN_VERSION}")
  endif()

  if ( NOT "${CMAKE_CFG_INTDIR}" STREQUAL "." )
    set(HAS_CFG_INTDIR TRUE CACHE INTERNAL "Generator uses intermediate configuration directory" )
    message( STATUS "Generator uses intermediate configuration directory: ${CMAKE_CFG_INTDIR}" )
  endif()

endif()
