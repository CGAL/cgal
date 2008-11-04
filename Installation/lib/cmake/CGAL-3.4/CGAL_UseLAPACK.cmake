# Do not check LAPACK_FOUND as it may be set by FindLAPACK.cmake
# if ( NOT LAPACK_FOUND )

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
    if ( "${LAPACK_DEFINITIONS}" MATCHES ".*LAPACK_USE_F2C.*" )
      add_definitions( "-DCGAL_USE_F2C" )
    endif()

    if ( BUILD_SHARED_LIBS )
      uniquely_add_flags( CMAKE_SHARED_LINKER_FLAGS ${LAPACK_LINKER_FLAGS} )
    else()
      uniquely_add_flags( CMAKE_MODULE_LINKER_FLAGS ${LAPACK_LINKER_FLAGS} )
    endif()

    link_directories( ${LAPACK_LIBRARIES_DIR} )
    link_libraries  ( ${LAPACK_LIBRARIES}     )

    # LAPACK requires BLAS
    include( CGAL_UseBLAS )

  endif()

# endif()

