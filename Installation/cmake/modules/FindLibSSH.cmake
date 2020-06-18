# - Try to find the LibSSH libraries
# This module defines:
#  LIBSSH_FOUND             - system has LibSSH lib
#  LIBSSH_INCLUDE_DIR       - the LibSSH include directory
#  LIBSSH_LIBRARIES_DIR     - directory where the LibSSH libraries are located
#  LIBSSH_LIBRARIES         - Link these to use LibSSH


include(FindPackageHandleStandardArgs)
include(${CMAKE_CURRENT_LIST_DIR}/CGAL_GeneratorSpecificSettings.cmake)

if(LIBSSH_INCLUDE_DIR)
  set(LIBSSH_in_cache TRUE)
else()
  set(LIBSSH_in_cache FALSE)
endif()
if(NOT LIBSSH_LIBRARIES)
  set(LIBSSH_in_cache FALSE)
endif()

# Is it already configured?
if( NOT LIBSSH_in_cache )

  find_path(LIBSSH_INCLUDE_DIR
    NAMES "libssh/libssh.h"
    )

  find_library(LIBSSH_LIBRARIES NAMES ssh libssh
    HINTS "/usr/lib"
    "usr/lib/x86_64-linux-gnu"
    PATH_SUFFIXES lib
    DOC "Path to the LIBSSH library"
    )
endif()

SET(LIBSSH_FOUND TRUE)
if( NOT LIBSSH_LIBRARIES OR NOT LIBSSH_INCLUDE_DIR)
  SET(LIBSSH_FOUND FALSE)
endif()
