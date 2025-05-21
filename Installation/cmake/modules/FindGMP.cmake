# - Try to find the GMP libraries
# This module defines:
#  GMP_FOUND             - system has GMP lib
#  GMP_INCLUDE_DIR       - the GMP include directory
#  GMP_LIBRARIES_DIR     - directory where the GMP libraries are located
#  GMP_LIBRARIES         - Link these to use GMP

include(FindPackageHandleStandardArgs)
include(${CMAKE_CURRENT_LIST_DIR}/CGAL_GeneratorSpecificSettings.cmake)

if(GMP_INCLUDE_DIR)
  set(GMP_in_cache TRUE)
else()
  set(GMP_in_cache FALSE)
endif()
if(NOT GMP_LIBRARIES)
  set(GMP_in_cache FALSE)
endif()

# Is it already configured?
if( NOT GMP_in_cache )

  find_path(GMP_INCLUDE_DIR
            NAMES gmp.h
            HINTS ENV GMP_INC_DIR
                  ENV GMP_DIR
                  $ENV{GMP_DIR}/include
                  ${CGAL_INSTALLATION_PACKAGE_DIR}/auxiliary/gmp/include
            PATH_SUFFIXES include
  	        DOC "The directory containing the GMP header files"
           )

  find_library(GMP_LIBRARY_RELEASE NAMES gmp libgmp-10 gmp-10 mpir
    HINTS ENV GMP_LIB_DIR
          ENV GMP_DIR
          $ENV{GMP_DIR}/lib
          ${CGAL_INSTALLATION_PACKAGE_DIR}/auxiliary/gmp/lib
    PATH_SUFFIXES lib
    DOC "Path to the Release GMP library"
    )

  find_library(GMP_LIBRARY_DEBUG NAMES gmpd gmp libgmp-10 gmp-10 mpir
    HINTS ENV GMP_LIB_DIR
          ENV GMP_DIR
          $ENV{GMP_DIR}/include
          ${CGAL_INSTALLATION_PACKAGE_DIR}/auxiliary/gmp/lib
    PATH_SUFFIXES lib
    DOC "Path to the Debug GMP library"
    )

  get_property(IS_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
  if(IS_MULTI_CONFIG)
    set(GMP_LIBRARIES debug ${GMP_LIBRARY_DEBUG} optimized ${GMP_LIBRARY_RELEASE})
  else()
    if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
      set(GMP_LIBRARIES ${GMP_LIBRARY_DEBUG})
    else()
      set(GMP_LIBRARIES ${GMP_LIBRARY_RELEASE})
    endif()
  endif()

  # Attempt to load a user-defined configuration for GMP if couldn't be found
  if ( NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARIES)
    include( GMPConfig OPTIONAL )
  endif()

endif()

find_package_handle_standard_args(GMP "DEFAULT_MSG" GMP_LIBRARIES GMP_INCLUDE_DIR)
