# This module setups the compiler for the RS3 library.
# It assumes that find_package(RS3) was already called.

if( RS3_FOUND AND NOT RS3_SETUP )

    include( CGAL_UseRS )

    # add rs3 parameters, if necessary (rs3 must be always after rsexport)
    message( STATUS "UseRS3" )
    message( STATUS "RS3 include:       ${RS3_INCLUDE_DIR}" )
    message( STATUS "RS3 definitions:   ${RS3_DEFINITIONS}" )
    message( STATUS "RS3 libraries:     ${RS3_LIBRARIES}" )

    include_directories ( SYSTEM ${RS3_INCLUDE_DIR} )
    add_definitions( ${RS3_DEFINITIONS} "-DCGAL_USE_RS3" )
    link_libraries( ${RS3_LIBRARIES} )

    # extract RS3 version from the file rsversion.h (based on Fernando
    # Cacciola's code for FindBoost.cmake)
    if( EXISTS "${RS3_INCLUDE_DIR}/rsversion.h" )
      FILE(READ "${RS3_INCLUDE_DIR}/rsversion.h" _rsversion_h_CONTENTS)
      STRING(REGEX REPLACE ".*#define[ \t]+rs_major[ \t]+([0-9]+).*"
                           "\\1" RS3_MAJOR "${_rsversion_h_CONTENTS}")
      STRING(REGEX REPLACE ".*#define[ \t]+rs_middle[ \t]+([0-9]+).*"
                           "\\1" RS3_MIDDLE "${_rsversion_h_CONTENTS}")
      STRING(REGEX REPLACE ".*#define[ \t]+rs_minor[ \t]+([0-9]+).*"
                           "\\1" RS3_MINOR "${_rsversion_h_CONTENTS}")
      SET( RS3_LIB_VERSION "${RS3_MAJOR}.${RS3_MIDDLE}.${RS3_MINOR}" )
      IS_VERSION_LESS( "${RS3_MAJOR}.${RS3_MIDDLE}" "3.1" RS3_OLD_INCLUDES )
    else( EXISTS "${RS3_INCLUDE_DIR}/rsversion.h" )
      # rsversion.h did not exist in old versions
      SET( RS3_LIB_VERSION "unknown" )
      SET( RS3_OLD_INCLUDES TRUE )
    endif( EXISTS "${RS3_INCLUDE_DIR}/rsversion.h" )

    message( STATUS "RS3 version is ${RS3_LIB_VERSION}" )

    if( RS3_OLD_INCLUDES )
      add_definitions( "-DCGAL_RS_OLD_INCLUDES" )
      message( STATUS "Using old RS signatures" )
    endif( RS3_OLD_INCLUDES )

    include(CGAL_UseMPFI)

  set (RS3_SETUP TRUE)

  endif( RS3_FOUND AND NOT RS3_SETUP )
