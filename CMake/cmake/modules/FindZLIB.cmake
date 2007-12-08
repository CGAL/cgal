# - Find ZLIB
# Find the native ZLIB includes and library
#
#  ZLIB_INCLUDE_DIR - where to find zlib.h, etc.
#  ZLIB_LIBRARIES   - List of libraries when using ZLIB.
#  ZLIB_FOUND       - True if ZLIB found.


IF (ZLIB_INCLUDE_DIR)
    # Already in cache, be silent
    SET(ZLIB_FIND_QUIETLY TRUE)
ENDIF (ZLIB_INCLUDE_DIR)

# After searching in standard places,
# search for precompiled ZLIB included with CGAL on Windows
SET(ZLIB_INCLUDE_DIR_SEARCH   /usr/local/include /usr/include)
SET(ZLIB_LIBRARIES_DIR_SEARCH /usr/lib /usr/local/lib)
IF(WIN32)
    SET(ZLIB_INCLUDE_DIR_SEARCH   ${ZLIB_INCLUDE_DIR_SEARCH}   ${CGAL_SOURCE_DIR}/auxiliary/zlib/include)
    SET(ZLIB_LIBRARIES_DIR_SEARCH ${ZLIB_LIBRARIES_DIR_SEARCH} ${CGAL_SOURCE_DIR}/auxiliary/zlib/lib)
ENDIF(WIN32)

FIND_PATH   (ZLIB_INCLUDE_DIR zlib.h 
       PATHS ${ZLIB_INCLUDE_DIR_SEARCH}
	     DOC "The directory containing the ZLIB include files")

SET(ZLIB_NAMES z zlib zdll)
FIND_LIBRARY(ZLIB_LIBRARY NAMES ${ZLIB_NAMES}
             PATHS ${ZLIB_LIBRARIES_DIR_SEARCH}
	           DOC "Path to the ZLIB library")

IF (ZLIB_INCLUDE_DIR AND ZLIB_LIBRARY)
    SET(ZLIB_FOUND TRUE)
    SET( ZLIB_LIBRARIES ${ZLIB_LIBRARY} )
ELSE (ZLIB_INCLUDE_DIR AND ZLIB_LIBRARY)
   SET(ZLIB_FOUND FALSE)
   SET( ZLIB_LIBRARIES )
ENDIF (ZLIB_INCLUDE_DIR AND ZLIB_LIBRARY)

# Print success/error message
if(ZLIB_FOUND)
    if(NOT ZLIB_FIND_QUIETLY)
        message(STATUS "Found ZLIB: ${ZLIB_LIBRARIES}")
    endif(NOT ZLIB_FIND_QUIETLY)
else(ZLIB_FOUND)
    IF(ZLIB_FIND_REQUIRED)
	MESSAGE(FATAL_ERROR "Could NOT find ZLIB. Set the ZLIB_INCLUDE_DIR and ZLIB_LIBRARIES cmake cache entries.")
    ELSE(ZLIB_FIND_REQUIRED)
	if(NOT ZLIB_FIND_QUIETLY)
	    MESSAGE(STATUS "Could NOT find ZLIB. Set the ZLIB_INCLUDE_DIR and ZLIB_LIBRARIES cmake cache entries.")
	endif(NOT ZLIB_FIND_QUIETLY)
    ENDIF(ZLIB_FIND_REQUIRED)
endif(ZLIB_FOUND)

MARK_AS_ADVANCED(ZLIB_LIBRARY ZLIB_INCLUDE_DIR)
