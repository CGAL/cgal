# LAPACK requires BLAS
include( CGAL_SetupBLAS )

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

    #get_dependency_version(LAPACK)

    add_definitions( ${LAPACK_DEFINITIONS} "-DCGAL_USE_LAPACK" )
    set( CGAL_3RD_PARTY_DEFINITIONS     ${CGAL_3RD_PARTY_DEFINITIONS} ${LAPACK_DEFINITIONS} )

    link_directories( ${LAPACK_LIBRARIES_DIR} )
    set( CGAL_3RD_PARTY_LIBRARIES_DIRS  ${CGAL_3RD_PARTY_LIBRARIES_DIRS} ${LAPACK_LIBRARIES_DIR} )
    set( CGAL_3RD_PARTY_LIBRARIES       ${CGAL_3RD_PARTY_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACK_LINKER_FLAGS} )

  endif()

endif()

