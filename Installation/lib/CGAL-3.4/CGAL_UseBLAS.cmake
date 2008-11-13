if ( BLAS_FOUND )

  message( STATUS "BLAS include:     ${BLAS_INCLUDE_DIR}" )
  include_directories ( ${BLAS_INCLUDE_DIR} )

  message( STATUS "BLAS definitions: ${BLAS_DEFINITIONS}" )
  add_definitions( ${BLAS_DEFINITIONS} "-DCGAL_USE_BLAS" )
  if ( "${BLAS_DEFINITIONS}" MATCHES ".*BLAS_USE_F2C.*" )
    add_definitions( "-DCGAL_USE_F2C" )
  endif()

  if (BLAS_LIBRARIES_DIR)
    message( STATUS "BLAS library directories:  ${BLAS_LIBRARIES_DIR}" )
    link_directories( ${BLAS_LIBRARIES_DIR} )
  endif()
  if (BLAS_LIBRARIES)
    message( STATUS "BLAS libraries:   ${BLAS_LIBRARIES}" )
    link_libraries  ( ${BLAS_LIBRARIES} )
  endif()

  message( STATUS "BLAS link flags:  ${BLAS_LINKER_FLAGS}" )
  if ( BUILD_SHARED_LIBS )
    uniquely_add_flags( CMAKE_SHARED_LINKER_FLAGS ${BLAS_LINKER_FLAGS} )
  else()
    uniquely_add_flags( CMAKE_MODULE_LINKER_FLAGS ${BLAS_LINKER_FLAGS} )
  endif()
  
endif()
