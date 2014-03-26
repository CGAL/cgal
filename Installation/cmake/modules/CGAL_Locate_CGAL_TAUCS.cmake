# Find TAUCS library shipped with CGAL
#
# This module searches for TAUCS in CGAL "auxiliary" folder
# and in in $CGAL_TAUCS_DIR environment variable.
#
# This module sets the following variables:
#  CGAL_TAUCS_FOUND - set to true if TAUCS library shipped with CGAL
#    is found
#  CGAL_TAUCS_PLATFORM - name of TAUCS subfolder corresponding to the current compiler
#  CGAL_TAUCS_INCLUDE_DIR - list of folders (using full path name) containing
#    TAUCS (and optionaly BLAS and LAPACK) headers
#  CGAL_TAUCS_LIBRARIES_DIR -list of folders (using full path name) containing
#    TAUCS (and optionaly BLAS and LAPACK) libraries

include(CGAL_Macros)

if ( NOT CGAL_TAUCS_FOUND )

  #
  # Find out TAUCS name for the current platform.
  # This code is a translation of TAUCS "configure" script.
  #

  # The first task is to figure out CMAKE_SYSTEM_NAME
  # (on unix this is uname -s, for windows it is Windows).
  #message("DEBUG: CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
  #message("DEBUG: CMAKE_SYSTEM_PROCESSOR = ${CMAKE_SYSTEM_PROCESSOR}")
  set( CGAL_TAUCS_PLATFORM "${CMAKE_SYSTEM_NAME}" )

  # Convert to lower case
  STRING(TOLOWER "${CGAL_TAUCS_PLATFORM}" CGAL_TAUCS_PLATFORM)

  # Sometimes uname returns a value that is
  # inconsistent with the way CGAL_TAUCS_PLATFORM is set. For example, on
  # Solaris, CGAL_TAUCS_PLATFORM=solaris but uname returns SunOS.
  if ( ${CGAL_TAUCS_PLATFORM} STREQUAL "sunos" )
    set( CGAL_TAUCS_PLATFORM "solaris" )
  endif()
  if ( ${CGAL_TAUCS_PLATFORM} STREQUAL "windows" )
    set( CGAL_TAUCS_PLATFORM "win32" )
  endif()

  # LS 2007: added "darwin_intel" for Intel Macs.
  # "darwin" = original Darwin platform = PowerPC architecture.
  if ( ${CGAL_TAUCS_PLATFORM} STREQUAL "darwin" )
    # CMAKE_SYSTEM_PROCESSOR=uname -p
    if ( ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "i386" )
      set( CGAL_TAUCS_PLATFORM "darwin_intel" )
    endif()
  endif()

  # LS 2007: append "64" if 64 bits processor (tested on Linux only)

  if ( ${CGAL_TAUCS_PLATFORM} STREQUAL "linux" )
    # CMAKE_SYSTEM_PROCESSOR=uname -p
    if ( ${CMAKE_SYSTEM_PROCESSOR} MATCHES ".*64.*" )
      set( CGAL_TAUCS_PLATFORM "${CGAL_TAUCS_PLATFORM}64" )
    endif()
  endif()

  #message("DEBUG: CGAL_TAUCS_PLATFORM = ${CGAL_TAUCS_PLATFORM}")


  #
  # Search for TAUCS folder.
  #

  #message("DEBUG: CGAL_INSTALLATION_PACKAGE_DIR = ${CGAL_INSTALLATION_PACKAGE_DIR}")

  # VC++ uses auto-link, thus we search for a folder containing the output of
  # build_taucs_win32_for_CGAL.bat/build_taucs_win64_for_CGAL.bat,
  # ie TAUCS libraries compiled for all Windows runtimes.
  if ( MSVC )

    # Check $CGAL_TAUCS_DIR environment variable
    fetch_env_var(CGAL_TAUCS_DIR)
    #message("DEBUG: CGAL_TAUCS_DIR = ${CGAL_TAUCS_DIR}")
    if (NOT "${CGAL_TAUCS_DIR}" STREQUAL "" AND EXISTS ${CGAL_TAUCS_DIR})

      set( CGAL_TAUCS_INCLUDE_DIR   "${CGAL_TAUCS_DIR}/include")
      set( CGAL_TAUCS_LIBRARIES_DIR "${CGAL_TAUCS_DIR}/lib"    )
      set( CGAL_TAUCS_FOUND TRUE )

    # Else, search for TAUCS in CGAL "auxiliary" folder
    elseif ( EXISTS "${CGAL_INSTALLATION_PACKAGE_DIR}/auxiliary/taucs/include/taucs.h" AND
             EXISTS "${CGAL_INSTALLATION_PACKAGE_DIR}/auxiliary/taucs/lib" )

      set( CGAL_TAUCS_INCLUDE_DIR   "${CGAL_INSTALLATION_PACKAGE_DIR}/auxiliary/taucs/include")
      set( CGAL_TAUCS_LIBRARIES_DIR "${CGAL_INSTALLATION_PACKAGE_DIR}/auxiliary/taucs/lib"    )
      set( CGAL_TAUCS_FOUND TRUE )

    endif()

  # On Unix, search for TAUCS "taucs_full" folder
  else ( MSVC )

    # Check $CGAL_TAUCS_DIR environment variable
    fetch_env_var(CGAL_TAUCS_DIR)
    #message("DEBUG: CGAL_TAUCS_DIR = ${CGAL_TAUCS_DIR}")
    if (NOT "${CGAL_TAUCS_DIR}" STREQUAL "" AND EXISTS ${CGAL_TAUCS_DIR})

        set( CGAL_TAUCS_INCLUDE_DIR   "${CGAL_TAUCS_DIR}/build/${CGAL_TAUCS_PLATFORM}"
                                      "${CGAL_TAUCS_DIR}/src" )
        set( CGAL_TAUCS_LIBRARIES_DIR "${CGAL_TAUCS_DIR}/external/lib/${CGAL_TAUCS_PLATFORM}"
                                      "${CGAL_TAUCS_DIR}/lib/${CGAL_TAUCS_PLATFORM}" )
        set( CGAL_TAUCS_FOUND TRUE )

    endif()

  endif ( MSVC )

  # Add variables to cache
  set( CGAL_TAUCS_DIR   "${CGAL_TAUCS_DIR}"
                        CACHE PATH "Directory containing BLAS, LAPACK and TAUCS" )

  #message("DEBUG: CGAL_TAUCS_INCLUDE_DIR = ${CGAL_TAUCS_INCLUDE_DIR}")
  #message("DEBUG: CGAL_TAUCS_LIBRARIES_DIR = ${CGAL_TAUCS_LIBRARIES_DIR}")
  #message("DEBUG: CGAL_TAUCS_FOUND = ${CGAL_TAUCS_FOUND}")

endif ( NOT CGAL_TAUCS_FOUND )

