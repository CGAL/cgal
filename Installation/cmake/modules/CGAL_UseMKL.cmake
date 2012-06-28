# This module setups the compiler for the MKL libraries.
# It assumes that find_package(MKL) was already called.

if ( MKL_FOUND AND NOT CGAL_MKL_SETUP )

  message( STATUS "MKL include:     ${MKL_INCLUDE_DIR}" )
  include_directories ( ${MKL_INCLUDE_DIR} )

  message( STATUS "MKL definitions: ${MKL_DEFINITIONS}" )
  add_definitions( ${MKL_DEFINITIONS} )
  if ( "${MKL_DEFINITIONS}" MATCHES ".*MKL_USE_F2C.*" )
    add_definitions( "-DCGAL_USE_F2C" )
  endif()

  if (MKL_LIBRARIES)
    message( STATUS "MKL libraries:   ${MKL_LIBRARIES}" )
    link_libraries( ${MKL_LIBRARIES} )
  endif()

  # Setup is done
  set ( CGAL_MKL_SETUP TRUE )

endif()

