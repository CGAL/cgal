if ( NOT BLAS_FOUND )

  find_package( BLAS )

  if ( BLAS_FOUND )

    message( STATUS "BLAS definitions: ${BLAS_DEFINITIONS}" )
    if (BLAS_LIBRARIES_DIR)
      message( STATUS "BLAS library directories:  ${BLAS_LIBRARIES_DIR}" )
    endif()
    if (BLAS_LIBRARIES)
      message( STATUS "BLAS libraries:   ${BLAS_LIBRARIES}" )
    endif()
    message( STATUS "BLAS link flags:  ${BLAS_LINKER_FLAGS}" )

    #get_dependency_version(BLAS)

    add_definitions( ${BLAS_DEFINITIONS} "-DCGAL_USE_BLAS" )
    set( CGAL_3RD_PARTY_DEFINITIONS     ${CGAL_3RD_PARTY_DEFINITIONS} ${BLAS_DEFINITIONS} )

    link_directories( ${BLAS_LIBRARIES_DIR} )
    set( CGAL_3RD_PARTY_LIBRARIES_DIRS  ${CGAL_3RD_PARTY_LIBRARIES_DIRS} ${BLAS_LIBRARIES_DIR} )
    set( CGAL_3RD_PARTY_LIBRARIES       ${CGAL_3RD_PARTY_LIBRARIES} ${BLAS_LIBRARIES} ${BLAS_LINKER_FLAGS} )

  endif()

endif()

