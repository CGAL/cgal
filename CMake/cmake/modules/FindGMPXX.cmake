# Try to find the GMPXX libraries
# GMPXX_FOUND - system has GMPXX lib
# GMPXX_INCLUDE_DIR - the GMPXX include directory
# GMPXX_LIBRARIES - Libraries needed to use GMPXX

# TODO: support Windows and MacOSX

# GMPXX needs GMP
find_package(GMP QUIET)
if(GMP_FOUND)
    if (GMPXX_INCLUDE_DIR AND GMPXX_LIBRARIES)
        # Already in cache, be silent
        set(GMPXX_FIND_QUIETLY TRUE)
    endif (GMPXX_INCLUDE_DIR AND GMPXX_LIBRARIES)

    find_path(GMPXX_INCLUDE_DIR NAMES gmpxx.h 
              PATHS ${GMP_INCLUDE_DIR_SEARCH}
      	      DOC "The directory containing the GMPXX include files")

    find_library(GMPXX_LIBRARIES NAMES gmpxx
                 PATHS ${GMP_LIBRARIES_DIR_SEARCH}
		             DOC "Path to the GMPXX library")

    if(GMPXX_INCLUDE_DIR AND GMPXX_LIBRARIES)
        set(GMPXX_FOUND TRUE)
    endif(GMPXX_INCLUDE_DIR AND GMPXX_LIBRARIES)

    # Print success/error message
    if(GMPXX_FOUND)
    	if(NOT GMPXX_FIND_QUIETLY)
    	    message(STATUS "Found GMPXX: ${GMPXX_LIBRARIES}")
    	endif(NOT GMPXX_FIND_QUIETLY)
        else(GMPXX_FOUND)
    	IF(GMPXX_FIND_REQUIRED)
    	    MESSAGE(FATAL_ERROR "Could NOT find GMPXX. Set the GMPXX_INCLUDE_DIR and GMPXX_LIBRARIES cmake cache entries.")
    	ELSE(GMPXX_FIND_REQUIRED)
    	    if(NOT GMPXX_FIND_QUIETLY)
    		MESSAGE(STATUS "Could NOT find GMPXX. Set the GMPXX_INCLUDE_DIR and GMPXX_LIBRARIES cmake cache entries.")
    	    endif(NOT GMPXX_FIND_QUIETLY)
    	ENDIF(GMPXX_FIND_REQUIRED)
    endif(GMPXX_FOUND)

    mark_as_advanced(GMPXX_INCLUDE_DIR GMPXX_LIBRARIES)
endif(GMP_FOUND)
