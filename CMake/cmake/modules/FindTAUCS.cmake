# Try to find the TAUCS libraries
# TAUCS_FOUND - system has TAUCS lib
# TAUCS_INCLUDE_DIR - the TAUCS include directory
# TAUCS_LIBRARIES - Libraries needed to use TAUCS

# TODO: support MacOSX

if (TAUCS_INCLUDE_DIR AND TAUCS_LIBRARIES)
  # Already in cache, be silent
  set(TAUCS_FIND_QUIETLY TRUE)
endif (TAUCS_INCLUDE_DIR AND TAUCS_LIBRARIES)

# After searching in standard places,
# search for precompiled TAUCS included with CGAL on Windows
if(WIN32)
    SET(TAUCS_INCLUDE_DIR_SEARCH   ${CGAL_SOURCE_DIR}/auxiliary/taucs/include)
    SET(TAUCS_LIBRARIES_DIR_SEARCH ${CGAL_SOURCE_DIR}/auxiliary/taucs/lib)
endif(WIN32)

find_path   (TAUCS_INCLUDE_DIR NAMES taucs.h
             PATHS ${TAUCS_INCLUDE_DIR_SEARCH}
      	     DOC "The directory containing the TAUCS include files")

if ( AUTO_LINK_ENABLED )
  if ( EXISTS "${TAUCS_LIBRARIES_DIR_SEARCH}" )
    SET(TAUCS_LIBRARIES ${TAUCS_LIBRARIES_DIR_SEARCH} )
  endif()  
else()
  find_library(TAUCS_LIBRARIES NAMES taucs PATHS
  	     ${TAUCS_LIBRARIES_DIR_SEARCH}
  	     DOC "Path to the TAUCS library")
endif()

if(TAUCS_INCLUDE_DIR AND TAUCS_LIBRARIES)
   set(TAUCS_FOUND TRUE)
endif()

# Print success/error message
if(TAUCS_FOUND)
    if(NOT TAUCS_FIND_QUIETLY)
        message(STATUS "Found TAUCS: ${TAUCS_LIBRARIES}")
    endif(NOT TAUCS_FIND_QUIETLY)
else(TAUCS_FOUND)
    if(TAUCS_FIND_REQUIRED)
	MESSAGE(FATAL_ERROR "Could NOT find TAUCS. Set the TAUCS_INCLUDE_DIR and TAUCS_LIBRARIES cmake cache entries.")
    else(TAUCS_FIND_REQUIRED)
	if(NOT TAUCS_FIND_QUIETLY)
	    MESSAGE(STATUS "Could NOT find TAUCS. Set the TAUCS_INCLUDE_DIR and TAUCS_LIBRARIES cmake cache entries.")
	endif(NOT TAUCS_FIND_QUIETLY)
    endif(TAUCS_FIND_REQUIRED)
endif(TAUCS_FOUND)

mark_as_advanced(TAUCS_INCLUDE_DIR TAUCS_LIBRARIES)
