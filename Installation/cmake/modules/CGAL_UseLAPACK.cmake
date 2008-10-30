# LAPACK requires BLAS
include( CGAL_UseBLAS )

if ( NOT LAPACK_FOUND )

  find_package( LAPACK )

  if ( LAPACK_FOUND )

    message( STATUS "LAPACK definitions: ${LAPACK_DEFINITIONS}" )
    if (LAPACK_LIBRARIES_DIR)
      message( STATUS "LAPACK library directories:  ${LAPACK_LIBRARIES_DIR}" )
    endif()
    if (LAPACK_LIBRARIES)
      message( STATUS "LAPACK libraries:   ${LAPACK_LIBRARIES}" )
    endif()
    message( STATUS "LAPACK link flags:  ${LAPACK_LINKER_FLAGS}" )

    add_definitions( ${LAPACK_DEFINITIONS} "-DCGAL_USE_LAPACK" )

    link_directories( ${LAPACK_LIBRARIES_DIR} )
    
    link_libraries( ${LAPACK_LIBRARIES} ${LAPACK_LINKER_FLAGS} )

  endif()

endif()

