# Try to find the CORE library shipped with CGAL
# CGAL_CORE_FOUND - system has CORE lib
# CGAL_CORE_INCLUDE_DIR - the CORE include directory
# CGAL_CORE_LIBRARIES - Libraries needed to use CORE

# TODO: support Windows and MacOSX

# CORE needs GMP
find_package(GMP QUIET)
if(GMP_FOUND)
    if (CGAL_CORE_INCLUDE_DIR AND CGAL_CORE_LIBRARIES)
        # Already in cache, be silent
        set(CGAL_CORE_FIND_QUIETLY TRUE)
    endif (CGAL_CORE_INCLUDE_DIR AND CGAL_CORE_LIBRARIES)

    # Find CORE include folder
    find_path(CGAL_CORE_INCLUDE_DIR NAMES CORE.h 
	      PATHS ${CGAL_SOURCE_DIR}/include/CORE
	      DOC "The directory containing the CORE include files shipped with CGAL")

    # We cannot search for the core++ library because it is not yet compiled
    # => hard code the name
    if (WIN32)
        set(CGAL_CORE_LIBRARIES ${CGAL_BINARY_DIR}/lib/core++.lib)
    else (WIN32)
        if(BUILD_SHARED_LIBS)
            set(CGAL_CORE_LIBRARIES ${CGAL_BINARY_DIR}/lib/libcore++.so)
        else(BUILD_SHARED_LIBS)
            set(CGAL_CORE_LIBRARIES ${CGAL_BINARY_DIR}/lib/libcore++.a)
        endif(BUILD_SHARED_LIBS)
    endif (WIN32)

    if(CGAL_CORE_INCLUDE_DIR AND CGAL_CORE_LIBRARIES)
        set(CGAL_CORE_FOUND TRUE)
    endif(CGAL_CORE_INCLUDE_DIR AND CGAL_CORE_LIBRARIES)

    # Print success/error message
    if(CGAL_CORE_FOUND)
	if(NOT CGAL_CORE_FIND_QUIETLY)
	    message(STATUS "Found CORE library shipped with CGAL: ${CGAL_CORE_LIBRARIES}")
	endif(NOT CGAL_CORE_FIND_QUIETLY)
    else(CGAL_CORE_FOUND)
	IF(CGAL_CORE_FIND_REQUIRED)
	    MESSAGE(FATAL_ERROR "Could NOT find CORE library shipped with CGAL. Set the CGAL_CORE_INCLUDE_DIR and CGAL_CORE_LIBRARIES cmake cache entries.")
	ELSE(CGAL_CORE_FIND_REQUIRED)
	    if(NOT CGAL_CORE_FIND_QUIETLY)
		MESSAGE(STATUS "Could NOT find CORE library shipped with CGAL. Set the CGAL_CORE_INCLUDE_DIR and CGAL_CORE_LIBRARIES cmake cache entries.")
	    endif(NOT CGAL_CORE_FIND_QUIETLY)
	ENDIF(CGAL_CORE_FIND_REQUIRED)
    endif(CGAL_CORE_FOUND)

    mark_as_advanced(CGAL_CORE_INCLUDE_DIR CGAL_CORE_LIBRARIES)
endif(GMP_FOUND)
