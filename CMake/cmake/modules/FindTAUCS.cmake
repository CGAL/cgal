# Try to find the TAUCS libraries
# TAUCS_FOUND - system has TAUCS lib
# TAUCS_INCLUDE_DIR - the TAUCS include directory
# TAUCS_LIBRARIES - Libraries needed to use TAUCS

# TODO: support MacOSX

if (TAUCS_INCLUDE_DIR AND TAUCS_LIBRARIES)
  # Already in cache, be silent
  set(TAUCS_FIND_QUIETLY TRUE)
endif()

SET(TAUCS_INCLUDE_DIR_SEARCH   ${CGAL_SOURCE_DIR}/auxiliary/taucs/include)
SET(TAUCS_LIBRARIES_DIR_SEARCH ${CGAL_SOURCE_DIR}/auxiliary/taucs/lib)

find_path(TAUCS_INCLUDE_DIR NAMES taucs.h
          PATHS ${TAUCS_INCLUDE_DIR_SEARCH}
      	  DOC "The directory containing the TAUCS include files"
         )

if ( AUTO_LINK_ENABLED )
  if ( EXISTS "${TAUCS_LIBRARIES_DIR_SEARCH}" )
    SET(TAUCS_LIBRARIES "")
  endif()  
else()
  find_library(TAUCS_LIBRARIES NAMES taucs PATHS
  	          ${TAUCS_LIBRARIES_DIR_SEARCH}
  	          DOC "Path to the TAUCS library"
              )
endif()

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(TAUCS "DEFAULT_MSG" TAUCS_INCLUDE_DIR TAUCS_LIBRARIES )
