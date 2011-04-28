# This module setups the compiler for the NTL library.
# It assumes that find_package(NTL) was already called.

if( NOT CGAL_NTL_SETUP )

  if( NTL_FOUND )

    message( STATUS "NTL include:        ${NTL_INCLUDE_DIR}" )
    message( STATUS "NTL definitions:    ${NTL_DEFINITIONS}" )
    message( STATUS "NTL libraries:      ${NTL_LIBRARIES}" )

    include_directories ( ${NTL_INCLUDE_DIRS} )
    add_definitions( ${NTL_DEFINITIONS} "-DCGAL_USE_NTL" )
    link_libraries( ${NTL_LIBRARIES} )

    set( CGAL_NTL_SETUP TRUE )

  endif( NTL_FOUND )


endif( NOT CGAL_NTL_SETUP )
