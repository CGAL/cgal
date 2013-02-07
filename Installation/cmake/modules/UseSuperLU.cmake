# This module setups the compiler for the SuperLU libraries.
# It assumes that find_package(SuperLU) was already called.

if ( SUPERLU_FOUND AND NOT CGAL_SUPERLU_SETUP )

  message( STATUS "SuperLU include:     ${SUPERLU_INCLUDES}" )
  include_directories ( ${SUPERLU_INCLUDES} )

  message( STATUS "SuperLU definitions: ${SUPERLU_DEFINITIONS}" )
  add_definitions( ${SUPERLU_DEFINITIONS} "-DCGAL_SUPERLU_ENABLED" )

  if (SUPERLU_LIBRARIES_DIR)
    message( STATUS "SuperLU library directories:  ${SUPERLU_LIBRARIES_DIR}" )
    link_directories( ${SUPERLU_LIBRARIES_DIR} )
  endif()
  if (SUPERLU_LIBRARIES)
    message( STATUS "SuperLU libraries:   ${SUPERLU_LIBRARIES}" )
    link_libraries( ${SUPERLU_LIBRARIES} )
  endif()

  # SuperLU requires BLAS and LAPACK
  include( ${LAPACK_USE_FILE} )

  # Setup is done
  set ( CGAL_SUPERLU_SETUP TRUE )

endif()

