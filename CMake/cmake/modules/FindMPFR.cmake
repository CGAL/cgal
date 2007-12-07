# Try to find the MPFR libraries
# MPFR_FOUND - system has MPFR lib
# MPFR_INCLUDE_DIR - the MPFR include directory
# MPFR_LIBRARIES - Libraries needed to use MPFR

# TODO: support Windows and MacOSX

# MPFR needs GMP
find_package(GMP QUIET)
if(GMP_FOUND)
    if (MPFR_INCLUDE_DIR AND MPFR_LIBRARIES)
        # Already in cache, be silent
        set(MPFR_FIND_QUIETLY TRUE)
    endif (MPFR_INCLUDE_DIR AND MPFR_LIBRARIES)

    find_path(MPFR_INCLUDE_DIR NAMES mpfr.h
	      DOC "The directory containing the MPFR include files")

    find_library(MPFR_LIBRARIES NAMES mpfr 
		 DOC "Path to the MPFR library")

    if(MPFR_INCLUDE_DIR AND MPFR_LIBRARIES)
        set(MPFR_FOUND TRUE)
    endif(MPFR_INCLUDE_DIR AND MPFR_LIBRARIES)

    # Print success/error message
    if(MPFR_FOUND)
	if(NOT MPFR_FIND_QUIETLY)
	    message(STATUS "Found MPFR: ${MPFR_LIBRARIES}")
	endif(NOT MPFR_FIND_QUIETLY)
    else(MPFR_FOUND)
	IF(MPFR_FIND_REQUIRED)
	    MESSAGE(FATAL_ERROR "Could NOT find MPFR. Set the MPFR_INCLUDE_DIR and MPFR_LIBRARIES cmake cache entries.")
	ELSE(MPFR_FIND_REQUIRED)
	    if(NOT MPFR_FIND_QUIETLY)
		MESSAGE(STATUS "Could NOT find MPFR. Set the MPFR_INCLUDE_DIR and MPFR_LIBRARIES cmake cache entries.")
	    endif(NOT MPFR_FIND_QUIETLY)
	ENDIF(MPFR_FIND_REQUIRED)
    endif(MPFR_FOUND)

    mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARIES)
endif(GMP_FOUND)
