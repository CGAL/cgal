if ( TAUCS_FOUND )

  message( STATUS "TAUCS include:     ${TAUCS_INCLUDE_DIR}" )
  include_directories ( ${TAUCS_INCLUDE_DIR} )

  message( STATUS "TAUCS definitions: ${TAUCS_DEFINITIONS}" )
  add_definitions( ${TAUCS_DEFINITIONS} )

  if (TAUCS_LIBRARIES_DIR)
    message( STATUS "TAUCS library directories:  ${TAUCS_LIBRARIES_DIR}" )
    link_directories( ${TAUCS_LIBRARIES_DIR} )
  endif()
  if (TAUCS_LIBRARIES)
    message( STATUS "TAUCS libraries:   ${TAUCS_LIBRARIES}" )
    link_libraries  ( ${TAUCS_LIBRARIES}     )
  endif()

  # TAUCS requires BLAS and LAPACK
  include(CGAL_UseLAPACK)

endif()
