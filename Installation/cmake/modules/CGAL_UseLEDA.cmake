# This module setups the compiler for the LEDA libraries.
# It assumes that find_package(LEDA) was already called.

if ( LEDA_FOUND AND NOT LEDA_SETUP )

  message( STATUS "LEDA include:     ${LEDA_INCLUDE_DIR}" )
  include_directories ( SYSTEM ${LEDA_INCLUDE_DIR} )

  message( STATUS "LEDA definitions: ${LEDA_DEFINITIONS}" )
  add_definitions( ${LEDA_DEFINITIONS} )
  if ( "${LEDA_DEFINITIONS}" MATCHES ".*LEDA_USE_F2C.*" )
    add_definitions( "-DCGAL_USE_F2C" )
  endif()

  if (LEDA_LIBRARIES_DIR)
    message( STATUS "LEDA library directories:  ${LEDA_LIBRARIES_DIR}" )
    link_directories( ${LEDA_LIBRARIES_DIR} )
  endif()
  if (LEDA_LIBRARIES)
    message( STATUS "LEDA libraries:   ${LEDA_LIBRARIES}" )
    link_libraries( ${LEDA_LIBRARIES} )
  endif()

  if (LEDA_CGAL_NO_STRICT_ALIASING)
    message( STATUS "${LEDA_CGAL_NO_STRICT_ALIASING}" )
  endif()
  if (LEDA_CGAL_LINK_X11)
    message( STATUS "${LEDA_CGAL_LINK_X11}" )
  endif()

  uniquely_add_flags( CMAKE_CXX_FLAGS ${LEDA_CXX_FLAGS} )
  uniquely_add_flags( CMAKE_EXE_LINKER_FLAGS ${LEDA_LINKER_FLAGS} )
  uniquely_add_flags( CMAKE_SHARED_LINKER_FLAGS ${LEDA_LINKER_FLAGS} )
  uniquely_add_flags( CMAKE_MODULE_LINKER_FLAGS ${LEDA_LINKER_FLAGS} )

 # Setup is done
  set ( LEDA_SETUP TRUE )

endif()
