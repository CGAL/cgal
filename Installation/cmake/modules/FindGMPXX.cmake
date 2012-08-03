# - Try to find the GMPXX libraries
# This module defines:
#   GMPXX_FOUND        - system has GMPXX lib
#   GMPXX_INCLUDE_DIR  - the GMPXX include directory
#   GMPXX_LIBRARIES    - Libraries needed to use GMPXX

# TODO: support Windows and MacOSX

# GMPXX needs GMP

find_package( GMP QUIET )

if(WITH_GMP AND GMP_FOUND)

  if (GMPXX_INCLUDE_DIR AND GMPXX_LIBRARIES)
    # Already in cache, be silent
    set(GMPXX_FIND_QUIETLY TRUE)
  endif()

  find_path(GMPXX_INCLUDE_DIR NAMES gmpxx.h 
            PATHS ${GMP_INCLUDE_DIR_SEARCH}
            DOC "The directory containing the GMPXX include files"
           )

  find_library(GMPXX_LIBRARIES NAMES gmpxx
               PATHS ${GMP_LIBRARIES_DIR_SEARCH}
               DOC "Path to the GMPXX library"
               )
               
  include(FindPackageHandleStandardArgs)
  
  find_package_handle_standard_args(GMPXX "DEFAULT_MSG" GMPXX_LIBRARIES GMPXX_INCLUDE_DIR )

else()

  message( FATAL_ERROR "GMPXX needs GMP")

endif()
