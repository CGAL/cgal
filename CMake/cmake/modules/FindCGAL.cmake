# Try to find the CGAL libraries
# CGAL_FOUND        - system has CGAL lib
# CGAL_DEFINITIONS  - the compiler settings to be compatible with CGAL
# CGAL_INCLUDE_DIRS - the CGAL include directories (including third-party libraries)
# CGAL_LIBRARIES    - libraries needed to use CGAL (including third-party libraries)

if (CGAL_INCLUDE_DIRS AND CGAL_LIBRARIES)
  # Already in cache, be silent
  set(CGAL_FIND_QUIETLY TRUE)
endif (CGAL_INCLUDE_DIRS AND CGAL_LIBRARIES)

# TODO: search CGAL in ${CGAL_INSTALL_DIR} and ${CGAL_BINARY_DIR}
#
#find_path   (CGAL_INCLUDE_DIRS NAMES CGAL/basic.h
#	     DOC "The directories containing include files for CGAL and third-party libraries")
#
#find_library(CGAL_LIBRARIES NAMES CGAL
#	     DOC "Path to CGAL and third-party libraries")

# TODO: include UseCGAL.cmake to define ${CGAL_DEFINITIONS}, ...?

if(CGAL_INCLUDE_DIRS AND CGAL_LIBRARIES)
   set(CGAL_FOUND TRUE)
endif(CGAL_INCLUDE_DIRS AND CGAL_LIBRARIES)

# Print success/error message
if(CGAL_FOUND)
    if(NOT CGAL_FIND_QUIETLY)
        message(STATUS "Found CGAL: ${CGAL_LIBRARIES}")
    endif(NOT CGAL_FIND_QUIETLY)
else(CGAL_FOUND)
    IF(CGAL_FIND_REQUIRED)
	MESSAGE(FATAL_ERROR "Could NOT find CGAL. Set the CGAL_INCLUDE_DIR and CGAL_LIBRARIES cmake cache entries.")
    ELSE(CGAL_FIND_REQUIRED)
	if(NOT CGAL_FIND_QUIETLY)
	    MESSAGE(STATUS "Could NOT find CGAL. Set the CGAL_INCLUDE_DIR and CGAL_LIBRARIES cmake cache entries.")
	endif(NOT CGAL_FIND_QUIETLY)
    ENDIF(CGAL_FIND_REQUIRED)
endif(CGAL_FOUND)

mark_as_advanced(CGAL_INCLUDE_DIRS CGAL_LIBRARIES)
