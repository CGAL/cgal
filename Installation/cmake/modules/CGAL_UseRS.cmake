# This module setups the compiler for the RS library.
# It assumes that find_package(RS) was already called.

if( NOT CGAL_RS_SETUP )

  #option( WITH_RS3 "Search for RS3" OFF )

  if( RS_FOUND )

    message( STATUS "RS include:        ${RS_INCLUDE_DIR}" )
    message( STATUS "RS definitions:    ${RS_DEFINITIONS}" )
    message( STATUS "RS libraries:      ${RS_LIBRARIES}" )

    if(CMAKE_OSX_ARCHITECTURES STREQUAL "ppc")
      message( STATUS "TLS is not supported on this architecture" )
      add_definitions( "-DCGAL_RS_NO_TLS" )
    endif(CMAKE_OSX_ARCHITECTURES STREQUAL "ppc")

    # add rs3 parameters, if necessary (rs3 must be always before rsexport)
    if(WITH_RS3 AND RS3_FOUND)

      message( STATUS "RS3 include:       ${RS3_INCLUDE_DIR}" )
      message( STATUS "RS3 definitions:   ${RS3_DEFINITIONS}" )
      message( STATUS "RS3 libraries:     ${RS3_LIBRARIES}" )

      if(CMAKE_OSX_ARCHITECTURES STREQUAL "ppc")
        message( STATUS "RS3 is not supported on this architecture" )
      else(CMAKE_OSX_ARCHITECTURES STREQUAL "ppc")
        include_directories ( ${RS3_INCLUDE_DIR} )
        add_definitions( ${RS_DEFINITIONS} "-DCGAL_USE_RS3" )
        link_libraries( ${RS3_LIBRARIES} )
      endif(CMAKE_OSX_ARCHITECTURES STREQUAL "ppc")
    endif()

    include_directories ( ${RS_INCLUDE_DIR} )
    add_definitions( ${RS_DEFINITIONS} "-DCGAL_USE_RS" )
    link_libraries( ${RS_LIBRARIES} )

  endif( RS_FOUND )

  set( CGAL_RS_SETUP TRUE )

endif( NOT CGAL_RS_SETUP )
