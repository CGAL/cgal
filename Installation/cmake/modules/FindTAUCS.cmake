# This module finds the TAUCS libraries.
#
# This module sets the following variables:
#  TAUCS_FOUND - Set to true if headers and libraries are found
#  TAUCS_INCLUDE_DIR - Directories containing the TAUCS header files
#  TAUCS_DEFINITIONS - Compilation options to use TAUCS
#  TAUCS_LIBRARIES_DIR - Directories containing the TAUCS libraries.
#     May be null if TAUCS_LIBRARIES contains libraries name using full path.
#  TAUCS_LIBRARIES - TAUCS libraries name.
#     May be null if the compiler supports auto-link (e.g. VC++).
#  TAUCS_USE_FILE - The name of the cmake module to include to compile
#     applications or libraries using TAUCS.

include(CGAL_GeneratorSpecificSettings)

# TAUCS requires LAPACK
if(TAUCS_FIND_QUIETLY OR NOT TAUCS_FIND_REQUIRED)
  find_package(LAPACK)
else()
  find_package(LAPACK REQUIRED)
endif()

if (NOT LAPACK_FOUND)

  message(STATUS "TAUCS requires LAPACK and BLAS.")
  set(TAUCS_FOUND FALSE)

# Is it already configured?
elseif (TAUCS_INCLUDE_DIR AND TAUCS_LIBRARIES_DIR)

  set(TAUCS_FOUND TRUE)

elseif (TAUCS_INCLUDE_DIR AND TAUCS_LIBRARIES)

  set(TAUCS_FOUND TRUE)

else()

  # reset variables
  set( TAUCS_DEFINITIONS "" ) # unused (yet)
  set( TAUCS_LIBRARIES "" )
  set( TAUCS_LIBRARIES_DIR "" )

  # Look first for the TAUCS library distributed with CGAL in auxiliary/taucs.
  # Set CGAL_TAUCS_FOUND, CGAL_TAUCS_INCLUDE_DIR and CGAL_TAUCS_LIBRARIES_DIR.
  include(CGAL_Locate_CGAL_TAUCS)

  # Search for TAUCS headers in ${CGAL_TAUCS_INCLUDE_DIR} (TAUCS shipped with CGAL),
  # else in $TAUCS_INC_DIR environment variable.
  if( CGAL_TAUCS_FOUND )
    set( TAUCS_INCLUDE_DIR    "${CGAL_TAUCS_INCLUDE_DIR}" )
  else()
    find_path(TAUCS_INCLUDE_DIR
              NAMES taucs.h
              HINTS ${CGAL_TAUCS_INCLUDE_DIR}
                    ENV TAUCS_DIR
                    ENV TAUCS_INC_DIR
              PATH_SUFFIXES taucs include
             )
  endif()

  # Search for TAUCS libraries in ${CGAL_TAUCS_LIBRARIES_DIR} (TAUCS shipped with CGAL),
  # else in $TAUCS_LIB_DIR environment variable.
  if( CGAL_TAUCS_FOUND AND CGAL_AUTO_LINK_ENABLED )
    # if VC++: done
    set( TAUCS_LIBRARIES_DIR  "${CGAL_TAUCS_LIBRARIES_DIR}" )
  else()
    find_library(TAUCS_LIBRARY
                 NAMES "taucs"
                 PATHS ${CGAL_TAUCS_LIBRARIES_DIR}
                       ENV TAUCS_DIR
                       ENV TAUCS_LIB_DIR
                 PATH_SUFFIXES taucs lib
                 DOC "TAUCS library"
                )
    find_library(METIS_LIBRARY
                 NAMES "metis"
                 PATHS ${CGAL_TAUCS_LIBRARIES_DIR}
                       ENV TAUCS
                       ENV TAUCS_LIB_DIR
                 PATH_SUFFIXES taucs lib
                 DOC "Metis library"
                )
    if(TAUCS_LIBRARY AND METIS_LIBRARY)
      set( TAUCS_LIBRARIES  "${TAUCS_LIBRARY};${METIS_LIBRARY}" )
    endif()
  endif()

  if (TAUCS_INCLUDE_DIR AND TAUCS_LIBRARIES_DIR)
    set(TAUCS_FOUND TRUE)
  elseif (TAUCS_INCLUDE_DIR AND TAUCS_LIBRARIES)
    set(TAUCS_FOUND TRUE)
  else()
    set(TAUCS_FOUND FALSE)
  endif()

  if(NOT TAUCS_FIND_QUIETLY)
    if(TAUCS_FOUND)
      message(STATUS "TAUCS libraries found.")
    else(TAUCS_FOUND)
      if(TAUCS_FIND_REQUIRED)
        message(FATAL_ERROR "TAUCS libraries not found. Please specify libraries location.")
      else()
        message(STATUS "TAUCS libraries not found. Please specify libraries location.")
      endif()
    endif(TAUCS_FOUND)
  endif(NOT TAUCS_FIND_QUIETLY)

  # Add variables to cache
  set( TAUCS_INCLUDE_DIR   "${TAUCS_INCLUDE_DIR}"
                            CACHE PATH "Directories containing the TAUCS header files" FORCE )
  set( TAUCS_DEFINITIONS   "${TAUCS_DEFINITIONS}"
                            CACHE STRING "Compilation options to use TAUCS" FORCE )
  set( TAUCS_LIBRARIES     "${TAUCS_LIBRARIES}"
                            CACHE FILEPATH "TAUCS libraries name" FORCE )
  set( TAUCS_LIBRARIES_DIR "${TAUCS_LIBRARIES_DIR}"
                            CACHE PATH "Directories containing the TAUCS libraries" FORCE )

  #message("DEBUG: TAUCS_INCLUDE_DIR = ${TAUCS_INCLUDE_DIR}")
  #message("DEBUG: TAUCS_DEFINITIONS = ${TAUCS_DEFINITIONS}")
  #message("DEBUG: TAUCS_LIBRARIES = ${TAUCS_LIBRARIES}")
  #message("DEBUG: TAUCS_LIBRARIES_DIR = ${TAUCS_LIBRARIES_DIR}")
  #message("DEBUG: TAUCS_FOUND = ${TAUCS_FOUND}")

endif(NOT LAPACK_FOUND)

if(TAUCS_FOUND)
  set(TAUCS_USE_FILE "CGAL_UseTAUCS")
endif(TAUCS_FOUND)
