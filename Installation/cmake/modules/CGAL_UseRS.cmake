# This module setups the compiler for the RS library.
# It assumes that find_package(RS) was already called.

if( RS_FOUND AND NOT RS_SETUP )

    message( STATUS "UseRS" )
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

    include_directories ( SYSTEM ${RS_INCLUDE_DIR} )
    add_definitions( ${RS_DEFINITIONS} "-DCGAL_USE_RS" )
    link_libraries( ${RS_LIBRARIES} )

  set (RS_SETUP TRUE)

endif( RS_FOUND AND NOT RS_SETUP )
