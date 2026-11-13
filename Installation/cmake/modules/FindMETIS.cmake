# -*- mode: cmake -*-

#
# METIS Find Module for MSTK
# Shamelessly stolen from Amanzi open source code https://software.lanl.gov/ascem/trac
#
# Usage:
#    Control the search through METIS_DIR or setting environment variable
#    METIS_ROOT to the METIS installation prefix.
#
#    This module does not search default paths!
#
#    Following variables are set:
#    METIS_FOUND            (BOOL)       Flag indicating if METIS was found
#    METIS_INCLUDE_DIR      (PATH)       Path to the METIS include file
#    METIS_INCLUDE_DIRS     (LIST)       List of all required include files
#    METIS_LIBRARY_DIR      (PATH)       Path to the METIS library
#    METIS_LIBRARY          (FILE)       METIS library
#    METIS_LIBRARIES        (LIST)       List of all required METIS libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# The name the caller used in find_package() (Metis, METIS, ...). CMake sets
# <ThatName>_FIND_QUIETLY / _FIND_REQUIRED, so read the flags under that name.
if(DEFINED CMAKE_FIND_PACKAGE_NAME)
    set(metis_pkg "${CMAKE_FIND_PACKAGE_NAME}")
else()
    set(metis_pkg METIS)
endif()
set(metis_quiet "${${metis_pkg}_FIND_QUIETLY}")

# Print a diagnostic only when find_package() was not called with QUIET
function(metis_msg msg)
    if(NOT metis_quiet)
        message(STATUS "${msg}")
    endif()
endfunction()

if ( METIS_LIBRARIES AND METIS_INCLUDE_DIRS )

    # Do nothing. Variables are set. No need to search again

else(METIS_LIBRARIES AND METIS_INCLUDE_DIRS)

    # Cache variables
    if(METIS_DIR)
        set(METIS_DIR "${METIS_DIR}" CACHE PATH "Path to search for METIS include and library files")
    endif()

    if(METIS_INCLUDE_DIR)
        set(METIS_INCLUDE_DIR "${METIS_INCLUDE_DIR}" CACHE PATH "Path to search for METIS include files")
    endif()

    if(METIS_LIBRARY_DIR)
        set(METIS_LIBRARY_DIR "${METIS_LIBRARY_DIR}" CACHE PATH "Path to search for METIS library files")
    endif()


    # Search for include files
    # Search order preference:
    #  (1) METIS_INCLUDE_DIR - check existence of path AND if the include files exist
    #  (2) METIS_DIR/<include>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(metis_inc_names "metis.h")
    if (METIS_INCLUDE_DIR)

        if (EXISTS "${METIS_INCLUDE_DIR}")

            find_path(metis_test_include_path
                      NAMES ${metis_inc_names}
                      HINTS ${METIS_INCLUDE_DIR}
                      NO_DEFAULT_PATH)
            if(NOT metis_test_include_path)
                metis_msg("Can not locate ${metis_inc_names} in ${METIS_INCLUDE_DIR}")
            endif()
            set(METIS_INCLUDE_DIR "${metis_test_include_path}")

        else()
            metis_msg("METIS_INCLUDE_DIR=${METIS_INCLUDE_DIR} does not exist")
            set(METIS_INCLUDE_DIR "METIS_INCLUDE_DIR-NOTFOUND")
        endif()

   else()

# Metis sometimes puts the include files in a subdir called Lib

        set(metis_inc_suffixes "include" "Lib")
        if(METIS_DIR)

            if (EXISTS "${METIS_DIR}" )

                find_path(METIS_INCLUDE_DIR
                          NAMES ${metis_inc_names}
                          HINTS ${METIS_DIR}
                          PATH_SUFFIXES ${metis_inc_suffixes}
                          NO_DEFAULT_PATH)

            else()
                 metis_msg("METIS_DIR=${METIS_DIR} does not exist")
                 set(METIS_INCLUDE_DIR "METIS_INCLUDE_DIR-NOTFOUND")
            endif()


        else()

            find_path(METIS_INCLUDE_DIR
                      NAMES ${metis_inc_names}
                      PATH_SUFFIXES ${metis_inc_suffixes})

        endif()

    endif()

    if ( NOT METIS_INCLUDE_DIR )
        metis_msg("Can not locate METIS include directory")
    endif()

    # Search for libraries
    # Search order preference:
    #  (1) METIS_LIBRARY_DIR - check existence of path AND if the library file exists
    #  (2) METIS_DIR/<lib,Lib>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(metis_lib_names "metis")
    if (METIS_LIBRARY_DIR)

        if (EXISTS "${METIS_LIBRARY_DIR}")

            find_library(METIS_LIBRARY
                         NAMES ${metis_lib_names}
                         HINTS ${METIS_LIBRARY_DIR}
                         NO_DEFAULT_PATH)
        else()
            metis_msg("METIS_LIBRARY_DIR=${METIS_LIBRARY_DIR} does not exist")
            set(METIS_LIBRARY "METIS_LIBRARY-NOTFOUND")
        endif()

    else()

        list(APPEND metis_lib_suffixes "lib" "Lib")
        if(METIS_DIR)

            if (EXISTS "${METIS_DIR}" )

                find_library(METIS_LIBRARY
                             NAMES ${metis_lib_names}
                             HINTS ${METIS_DIR}
                             PATH_SUFFIXES ${metis_lib_suffixes}
                             NO_DEFAULT_PATH)

            else()
                 metis_msg("METIS_DIR=${METIS_DIR} does not exist")
                 set(METIS_LIBRARY "METIS_LIBRARY-NOTFOUND")
            endif()


        else()

            find_library(METIS_LIBRARY
                         NAMES ${metis_lib_names}
                         PATH_SUFFIXES ${metis_lib_suffixes})

        endif()

    endif()

    if ( NOT METIS_LIBRARY )
        metis_msg("Can not locate METIS library")
    endif()


    # Define prerequisite packages
    set(METIS_INCLUDE_DIRS ${METIS_INCLUDE_DIR})
    set(METIS_LIBRARIES    ${METIS_LIBRARY})


endif(METIS_LIBRARIES AND METIS_INCLUDE_DIRS )

# Use the caller's package name so QUIET/REQUIRED are honored. This sets both
# <Name>_FOUND (e.g. Metis_FOUND) and the upper-case METIS_FOUND.
find_package_handle_standard_args(${metis_pkg}
                                  REQUIRED_VARS METIS_LIBRARIES METIS_INCLUDE_DIRS)

mark_as_advanced(
  METIS_INCLUDE_DIR
  METIS_INCLUDE_DIRS
  METIS_LIBRARY
  METIS_LIBRARIES
  METIS_LIBRARY_DIR
)

unset(metis_pkg)
unset(metis_quiet)
