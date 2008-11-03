# Do not check TAUCS_FOUND as it may be set by FindTAUCS.cmake
# if ( NOT TAUCS_FOUND )

  find_package( TAUCS )

  if ( TAUCS_FOUND )

    message( STATUS "TAUCS include:     ${TAUCS_INCLUDE_DIR}" )
    message( STATUS "TAUCS definitions: ${TAUCS_DEFINITIONS}" )
    if (TAUCS_LIBRARIES_DIR)
      message( STATUS "TAUCS library directories:  ${TAUCS_LIBRARIES_DIR}" )
    endif()
    if (TAUCS_LIBRARIES)
      message( STATUS "TAUCS libraries:   ${TAUCS_LIBRARIES}" )
    endif()

    include_directories ( ${TAUCS_INCLUDE_DIR} )

    add_definitions( ${TAUCS_DEFINITIONS} "-DCGAL_USE_TAUCS" )

    link_directories( ${TAUCS_LIBRARIES_DIR} )
    link_libraries  ( ${TAUCS_LIBRARIES}     )

    # TAUCS requires BLAS and LAPACK
    include(CGAL_UseLAPACK)

  endif()

# endif()

