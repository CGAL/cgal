if (LEDA_INCLUDE_DIR AND LEDA_LIBRARIES )

  set(LEDA_FOUND TRUE)

else()

  find_path(LEDA_INCLUDE_DIR
            NAMES "LEDA/basic.h" "LEDA/system/basic.h"
            HINTS ENV LEDA_INC_DIR
                  ENV LEDA_DIR
            PATH_SUFFIXES incl
            DOC "The directory containing the LEDA header files WITHOUT the LEDA prefix"
          )

  find_library(LEDA_LIBRARY_RELEASE NAMES "leda_md" "leda" "leda_numbers_md" "leda_numbers"
               HINTS ENV LEDA_LIB_DIR
                     ENV LEDA_DIR
               PATH_SUFFIXES lib
               DOC "Path to the LEDA library"
              )

  find_library(LEDA_LIBRARY_DEBUG NAMES "leda_mdd" "ledaD" "leda_numbers_mdd" "leda_numbers"
               HINTS ENV LEDA_LIB_DIR
                     ENV LEDA_DIR
               PATH_SUFFIXES lib
               DOC "Path to the LEDA library"
              )

  if ( NOT LEDA_INCLUDE_DIR )
    typed_cache_set( PATH "The directory containing the LEDA header files WITHOUT the LEDA prefix" LEDA_INCLUDE_DIR "$ENV{LEDA_INC_DIR}" )
  endif()

  if ( NOT LEDA_DEFINITIONS )
    typed_cache_set( STRING "Definitions for the LEDA library" LEDA_DEFINITIONS "$ENV{LEDA_DEFINITIONS}" )
  endif()

  if ( NOT LEDA_CXX_FLAGS )
    typed_cache_set( STRING "Compiler flags for the LEDA library" LEDA_CXX_FLAGS "$ENV{LEDA_CXX_FLAGS}" )
  endif()

  if ( NOT LEDA_LINKER_FLAGS )
    typed_cache_set( STRING "Linker flags for the LEDA library" LEDA_LINKER_FLAGS "$ENV{LEDA_LINKER_FLAGS}" )
  endif()

  if ( NOT LEDA_LIBRARY_RELEASE )
    typed_cache_set( FILEPATH "The LEDA release-mode libraries" LEDA_LIBRARY_RELEASE "$ENV{LEDA_LIBRARY_RELEASE}" )
  endif()

  if ( NOT LEDA_LIBRARY_DEBUG )
    typed_cache_set( FILEPATH "The LEDA debug-mode libraries" LEDA_LIBRARY_DEBUG "$ENV{LEDA_LIBRARY_DEBUG}" )
  endif()

  if(LEDA_LIBRARY_RELEASE)
    if(LEDA_LIBRARY_DEBUG)
      set(LEDA_LIBRARIES_ optimized ${LEDA_LIBRARY_RELEASE} debug ${LEDA_LIBRARY_DEBUG})
    else()
      set(LEDA_LIBRARIES_ ${LEDA_LIBRARY_RELEASE})
    endif()
  endif()

  set(LEDA_LIBRARIES ${LEDA_LIBRARIES_} CACHE FILEPATH "The LEDA library")

endif()

set( LEDA_BASIC_H "${LEDA_INCLUDE_DIR}/LEDA/system/basic.h" )
if ( NOT EXISTS ${LEDA_BASIC_H} )
  set( LEDA_BASIC_H "${LEDA_INCLUDE_DIR}/LEDA/basic.h" )
endif()

if ( EXISTS ${LEDA_BASIC_H} )

  file(READ "${LEDA_BASIC_H}" LEDA_BASIC_H_CONTENTS)

  string(REGEX REPLACE ".*#define __LEDA__ ([0-9]+).*" "\\1" LEDA_VERSION "${LEDA_BASIC_H_CONTENTS}")

  message( STATUS "USING LEDA_VERSION = '${LEDA_VERSION}'" )

endif()

if ( NOT LEDA_VERSION AND NOT "$ENV{LEDA_VERSION}" STREQUAL "" )
  typed_cache_set( STRING "The version of the LEDA library" LEDA_VERSION "$ENV{LEDA_VERSION}" )
endif()

if ( LEDA_VERSION )
  if ( NOT "${LEDA_DEFINITIONS}" MATCHES "-DCGAL_LEDA_VERSION=${LEDA_VERSION}" )
    typed_cache_set( STRING "Definitions for the LEDA library" LEDA_DEFINITIONS "${LEDA_DEFINITIONS}" "-DCGAL_LEDA_VERSION=${LEDA_VERSION}" )
  endif()
endif()

if ( LEDA_INCLUDE_DIR AND LEDA_LIBRARIES)
  set(LEDA_FOUND TRUE)

  if ( CMAKE_COMPILER_IS_GNUCXX )
    get_dependency_version (GCC)
    if ( NOT "${GCC_VERSION}" VERSION_LESS "4.4" )
      set(LEDA_CGAL_NO_STRICT_ALIASING TRUE)
      typed_cache_set( INTERNAL "Add -fno-strict-aliasing on gcc >= 4.4" LEDA_CGAL_NO_STRICT_ALIASING "Using LEDA with gcc version 4.4 or later: Adding -fno-strict-aliasing")
      uniquely_add_flags (LEDA_CXX_FLAGS "-fno-strict-aliasing")
    endif()
  endif()
  if ( UNIX AND IS_DIRECTORY "${LEDA_INCLUDE_DIR}/LEDA/graphics")
    set(LEDA_CGAL_LINK_X11 TRUE)
    typed_cache_set( INTERNAL "Link against X11 on *nix" LEDA_CGAL_LINK_X11 "Using LEDA with gcc on *nix: Adding -lX11")
    uniquely_add_flags( LEDA_LINKER_FLAGS "-lX11" )
  endif()

  set(LEDA_USE_FILE "CGAL_UseLEDA")
endif()

