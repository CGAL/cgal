# Try to find the CORE libraries
# CORE_FOUND - system has CORE lib
# CORE_INCLUDE_DIR - the CORE include directory
# CORE_LIBRARIES - Libraries needed to use CORE

# TODO: support Windows and MacOSX

# CORE needs GMP
find_package(GMP QUIET)
if(GMP_FOUND)
    if (CORE_INCLUDE_DIR AND CORE_LIBRARIES)
        # Already in cache, be silent
        set(CORE_FIND_QUIETLY TRUE)
    endif (CORE_INCLUDE_DIR AND CORE_LIBRARIES)

    find_path(CORE_INCLUDE_DIR NAMES CORE.h
	      DOC "The directory containing the CORE include files")

    find_library(CORE_LIBRARIES NAMES core++
		 DOC "Path to the core++ library")

    if(CORE_INCLUDE_DIR AND CORE_LIBRARIES)
        set(CORE_FOUND TRUE)
    endif(CORE_INCLUDE_DIR AND CORE_LIBRARIES)

    # Print success/error message
    if(CORE_FOUND)
	if(NOT CORE_FIND_QUIETLY)
	    message(STATUS "Found CORE: ${CORE_LIBRARIES}")
	endif(NOT CORE_FIND_QUIETLY)
    else(CORE_FOUND)
	IF(CORE_FIND_REQUIRED)
	    MESSAGE(FATAL_ERROR "Could NOT find CORE. Set the CORE_INCLUDE_DIR and CORE_LIBRARIES cmake cache entries.")
	ELSE(CORE_FIND_REQUIRED)
	    if(NOT CORE_FIND_QUIETLY)
		MESSAGE(STATUS "Could NOT find CORE. Set the CORE_INCLUDE_DIR and CORE_LIBRARIES cmake cache entries.")
	    endif(NOT CORE_FIND_QUIETLY)
	ENDIF(CORE_FIND_REQUIRED)
    endif(CORE_FOUND)

    mark_as_advanced(CORE_INCLUDE_DIR CORE_LIBRARIES)
endif(GMP_FOUND)
