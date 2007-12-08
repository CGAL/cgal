# Try to find the GMP libraries
# GMP_FOUND - system has GMP lib
# GMP_INCLUDE_DIR - the GMP include directory
# GMP_LIBRARIES - Libraries needed to use GMP

# TODO: support MacOSX

if (GMP_INCLUDE_DIR AND GMP_LIBRARIES)
  # Already in cache, be silent
  set(GMP_FIND_QUIETLY TRUE)
endif (GMP_INCLUDE_DIR AND GMP_LIBRARIES)

# After searching in standard places,
# search for precompiled GMP included with CGAL on Windows
if(WIN32)
    SET(GMP_INCLUDE_DIR_SEARCH   ${CGAL_SOURCE_DIR}/auxiliary/gmp/include)
    SET(GMP_LIBRARIES_DIR_SEARCH ${CGAL_SOURCE_DIR}/auxiliary/gmp/lib)
endif(WIN32)

find_path   (GMP_INCLUDE_DIR NAMES gmp.h PATHS
	     ${GMP_INCLUDE_DIR_SEARCH}
	     DOC "The directory containing the GMP include files")

if ( AUTO_LINK_ENABLED )
  if ( EXISTS "${GMP_LIBRARIES_DIR_SEARCH}" )
    SET(GMP_LIBRARIES ${GMP_LIBRARIES_DIR_SEARCH} )
  endif()  
else()
  find_library(GMP_LIBRARIES NAMES gmp PATHS
  	     ${GMP_LIBRARIES_DIR_SEARCH}
  	     DOC "Path to the GMP library")
endif()

if(GMP_INCLUDE_DIR AND GMP_LIBRARIES)
   set(GMP_FOUND TRUE)
endif()

# Print success/error message
if(GMP_FOUND)
    if(NOT GMP_FIND_QUIETLY)
        message(STATUS "Found GMP: ${GMP_LIBRARIES}")
    endif(NOT GMP_FIND_QUIETLY)
else(GMP_FOUND)
    if(GMP_FIND_REQUIRED)
	MESSAGE(FATAL_ERROR "Could NOT find GMP. Set the GMP_INCLUDE_DIR and GMP_LIBRARIES cmake cache entries.")
    else(GMP_FIND_REQUIRED)
	if(NOT GMP_FIND_QUIETLY)
	    MESSAGE(STATUS "Could NOT find GMP. Set the GMP_INCLUDE_DIR and GMP_LIBRARIES cmake cache entries.")
	endif(NOT GMP_FIND_QUIETLY)
    endif(GMP_FIND_REQUIRED)
endif(GMP_FOUND)

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARIES)
