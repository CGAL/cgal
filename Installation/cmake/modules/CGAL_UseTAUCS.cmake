# TAUCS requires BLAS and LAPACK
include(CGAL_UseBLAS)
include(CGAL_UseLAPACK)

if ( NOT TAUCS_FOUND )

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
    lini_libraries  ( ${TAUCS_LIBRARIES}     )

  endif()

endif()

