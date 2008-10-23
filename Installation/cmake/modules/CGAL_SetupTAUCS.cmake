# TAUCS requires BLAS and LAPACK
include(CGAL_SetupBLAS)
include(CGAL_SetupLAPACK)

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

    #get_dependency_version(TAUCS)

    include_directories ( ${TAUCS_INCLUDE_DIR} )

    add_definitions( ${TAUCS_DEFINITIONS} "-DCGAL_USE_TAUCS" )
    set( CGAL_3RD_PARTY_DEFINITIONS     ${CGAL_3RD_PARTY_DEFINITIONS} ${TAUCS_DEFINITIONS} )

    link_directories( ${TAUCS_LIBRARIES_DIR} )
    set( CGAL_3RD_PARTY_LIBRARIES_DIRS  ${CGAL_3RD_PARTY_LIBRARIES_DIRS} ${TAUCS_LIBRARIES_DIR} )
    set( CGAL_3RD_PARTY_LIBRARIES       ${CGAL_3RD_PARTY_LIBRARIES} ${TAUCS_LIBRARIES} )

  endif()

endif()

