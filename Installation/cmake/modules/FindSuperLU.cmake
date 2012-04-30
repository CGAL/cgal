
# Umfpack lib usually requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.

find_package(LAPACK)

if (NOT LAPACK_FOUND)
  message(STATUS "SuperLU requires LAPACK and BLAS that could not be found.")
  set(SUPERLU_FOUND FALSE)
else()
  if (SUPERLU_INCLUDES AND SUPERLU_LIBRARIES)
    set(SUPERLU_FIND_QUIETLY TRUE)
  endif (SUPERLU_INCLUDES AND SUPERLU_LIBRARIES)

  find_path(SUPERLU_INCLUDES
    NAMES
    supermatrix.h
    HINTS
    $ENV{SUPERLUDIR}/SRC
    $ENV{SUPERLU_INC_DIR}
    PATHS
    ${INCLUDE_INSTALL_DIR}
    PATH_SUFFIXES
    superlu
    SRC
  )

  find_library(SUPERLU_LIBRARIES superlu
               HINTS
               $ENV{SUPERLU_LIB_DIR} 
               $ENV{SUPERLUDIR}/lib
               PATHS 
               ${LIB_INSTALL_DIR}
               PATH_SUFFIXES
               lib
  )
    
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(SUPERLU DEFAULT_MSG
                                    SUPERLU_INCLUDES SUPERLU_LIBRARIES)

  mark_as_advanced(SUPERLU_INCLUDES SUPERLU_LIBRARIES)
endif()

if(SUPERLU_FOUND)
  set(SUPERLU_USE_FILE "UseSuperLU")
endif(SUPERLU_FOUND)