# This module setups the compiler for the BLAS libraries.
# It assumes that find_package(BLAS) was already called.

if ( BLAS_FOUND AND NOT BLAS_SETUP )

  message( STATUS "UseBLAS" )
  message( STATUS "BLAS include:     ${BLAS_INCLUDE_DIR}" )
  include_directories ( SYSTEM ${BLAS_INCLUDE_DIR} )

  message( STATUS "BLAS definitions: ${BLAS_DEFINITIONS}" )
  add_definitions( ${BLAS_DEFINITIONS} )
  if ( "${BLAS_DEFINITIONS}" MATCHES ".*BLAS_USE_F2C.*" )
    add_definitions( "-DCGAL_USE_F2C" )
  endif()

  if (BLAS_LIBRARIES_DIR)
    message( STATUS "BLAS library directories:  ${BLAS_LIBRARIES_DIR}" )
    link_directories( ${BLAS_LIBRARIES_DIR} )
  endif()
  if (BLAS_LIBRARIES)
    message( STATUS "BLAS libraries:   ${BLAS_LIBRARIES}" )
    link_libraries( ${BLAS_LIBRARIES} )
  endif()

  message( STATUS "BLAS link flags:  ${BLAS_LINKER_FLAGS}" )
  if ( BUILD_SHARED_LIBS )
    uniquely_add_flags( CMAKE_SHARED_LINKER_FLAGS ${BLAS_LINKER_FLAGS} )
  else()
    uniquely_add_flags( CMAKE_MODULE_LINKER_FLAGS ${BLAS_LINKER_FLAGS} )
  endif()

  # Setup is done
  set ( BLAS_SETUP TRUE )

endif()

