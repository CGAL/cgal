# This module setups the compiler for the TAUCS libraries.
# It assumes that find_package(TAUCS) was already called.

if ( TAUCS_FOUND AND NOT TAUCS_SETUP )

  message( STATUS "UseTAUCS" )
  message( STATUS "TAUCS include:     ${TAUCS_INCLUDE_DIR}" )
  include_directories ( SYSTEM ${TAUCS_INCLUDE_DIR} )

  message( STATUS "TAUCS definitions: ${TAUCS_DEFINITIONS}" )
  add_definitions( ${TAUCS_DEFINITIONS} "-DCGAL_USE_TAUCS" )

  if (TAUCS_LIBRARIES_DIR)
    message( STATUS "TAUCS library directories:  ${TAUCS_LIBRARIES_DIR}" )
    link_directories( ${TAUCS_LIBRARIES_DIR} )
  endif()
  if (TAUCS_LIBRARIES)
    message( STATUS "TAUCS libraries:   ${TAUCS_LIBRARIES}" )
    link_libraries( ${TAUCS_LIBRARIES} )
  endif()

  # TAUCS requires BLAS and LAPACK
  include( ${LAPACK_USE_FILE} )

  # Setup is done
  set ( TAUCS_SETUP TRUE )

  add_definitions(-DCGAL_TAUCS_ENABLED)

endif()

