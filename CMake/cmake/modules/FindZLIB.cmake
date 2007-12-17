# - Find ZLIB
# Find the native ZLIB includes and library
#
#  ZLIB_FOUND       - True if ZLIB found.
#  ZLIB_INCLUDE_DIR - where to find zlib.h, etc.
#  ZLIB_LIBRARIES   - List of libraries when using ZLIB.

IF (ZLIB_INCLUDE_DIR AND ZLIB_LIBRARIES )
  # Already in cache, be silent
  SET(ZLIB_FIND_QUIETLY TRUE)
ENDIF()

SET(ZLIB_INCLUDE_DIR_SEARCH   ${CGAL_SOURCE_DIR}/auxiliary/zlib/include)
SET(ZLIB_LIBRARIES_DIR_SEARCH ${CGAL_SOURCE_DIR}/auxiliary/zlib/lib)

FIND_PATH(ZLIB_INCLUDE_DIR zlib.h 
          PATHS ${ZLIB_INCLUDE_DIR_SEARCH}
	        DOC "The directory containing the ZLIB include files"
         )

SET(ZLIB_NAMES z zlib zdll)
FIND_LIBRARY(ZLIB_LIBRARIES NAMES ${ZLIB_NAMES}
             PATHS ${ZLIB_LIBRARIES_DIR_SEARCH}
	           DOC "Path to the ZLIB library"
            )

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ZLIB "DEFAULT_MSG" ZLIB_LIBRARIES ZLIB_INCLUDE_DIR )
