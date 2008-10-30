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

    add_definitions( ${BLAS_DEFINITIONS} "-DCGAL_USE_BLAS" )

    link_directories( ${BLAS_LIBRARIES_DIR} )
    
    link_libraries( ${BLAS_LIBRARIES} ${BLAS_LINKER_FLAGS} )

  endif()

endif()

