# This module setups the compiler for the RS library.
# It assumes that find_package(RS) was already called.

if( NOT CGAL_RS_SETUP )

  if( RS_FOUND )

    message( STATUS "RS include:        ${RS_INCLUDE_DIR}" )
    message( STATUS "RS definitions:    ${RS_DEFINITIONS}" )
    message( STATUS "RS libraries:      ${RS_LIBRARIES}" )

    if( APPLE AND CMAKE_COMPILER_IS_GNUCXX )
      include( CGAL_VersionUtils )
      EXEC_PROGRAM( ${CMAKE_CXX_COMPILER}
                    ARGS -dumpversion
                    OUTPUT_VARIABLE RS_GXX_VERSION )
      VERSION_DECOMPOSE( ${RS_GXX_VERSION} GXX_MAJ GXX_MIN GXX_PAT GXX_TWE )
      IS_VERSION_LESS( "${GXX_MAJ}.${GXX_MIN}" "4.3" IS_OLD_GXX )
      if( IS_OLD_GXX )
        message( STATUS "TLS is not supported by g++<4.3 on Mac OS X" )
        add_definitions( "-DCGAL_RS_NO_TLS" )
      endif( IS_OLD_GXX )
    endif( APPLE AND CMAKE_COMPILER_IS_GNUCXX )

    include_directories ( ${RS_INCLUDE_DIR} )
    add_definitions( ${RS_DEFINITIONS} "-DCGAL_USE_RS" )
    link_libraries( ${RS_LIBRARIES} )

    # add rs3 parameters, if necessary (rs3 must be always after rsexport)
    if( RS3_FOUND )

      message( STATUS "RS3 include:       ${RS3_INCLUDE_DIR}" )
      message( STATUS "RS3 definitions:   ${RS3_DEFINITIONS}" )
      message( STATUS "RS3 libraries:     ${RS3_LIBRARIES}" )

      include_directories ( ${RS3_INCLUDE_DIR} )
      add_definitions( ${RS_DEFINITIONS} "-DCGAL_USE_RS3" )
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
        IS_VERSION_LESS( "${RS3_MAJOR}.${RS3_MIDDLE}" "3.1" RS_OLD_INCLUDES )
      else( EXISTS "${RS3_INCLUDE_DIR}/rsversion.h" )
        # rsversion.h did not exist in old versions
        SET( RS3_LIB_VERSION "unknown" )
        SET( RS_OLD_INCLUDES TRUE )
      endif( EXISTS "${RS3_INCLUDE_DIR}/rsversion.h" )

      message( STATUS "RS version is ${RS3_LIB_VERSION}" )

      if( RS_OLD_INCLUDES )
        add_definitions( "-DCGAL_RS_OLD_INCLUDES" )
        message( STATUS "Using old RS signatures" )
      endif( RS_OLD_INCLUDES )

    endif( RS3_FOUND )

  endif( RS_FOUND )

  set( CGAL_RS_SETUP TRUE )

endif( NOT CGAL_RS_SETUP )
