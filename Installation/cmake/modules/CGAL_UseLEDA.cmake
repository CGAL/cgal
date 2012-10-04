# This module setups the compiler for the LEDA libraries.
# It assumes that find_package(LEDA) was already called.

if ( LEDA_FOUND AND NOT LEDA_SETUP )

  message( STATUS "UseLEDA" )
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

  if ( CMAKE_COMPILER_IS_GNUCXX )
    get_dependency_version (GCC)
    if ( NOT "${GCC_VERSION}" VERSION_LESS "4.1" )
      message( STATUS "Using LEDA with gcc version 4.1 or later. Adding -ffriend-injection" )
      uniquely_add_flags (CMAKE_CXX_FLAGS "-ffriend-injection")
    endif()
    if ( NOT "${GCC_VERSION}" VERSION_LESS "4.4" )
      message( STATUS "Using LEDA with gcc version 4.4 or later. Adding -fno-strict-aliasing" )
      uniquely_add_flags (CMAKE_CXX_FLAGS "-fno-strict-aliasing")
    endif()
    if ( UNIX )
      message( STATUS "Using LEDA with gcc on *nix. Adding -lX11" )
      uniquely_add_flags( CMAKE_SHARED_LINKER_FLAGS "-lX11" )
      uniquely_add_flags( CMAKE_MODULE_LINKER_FLAGS "-lX11" )
    endif()
  endif()

  uniquely_add_flags( CMAKE_CXX_FLAGS ${LEDA_CXX_FLAGS} )
  uniquely_add_flags( CMAKE_SHARED_LINKER_FLAGS ${LEDA_LINKER_FLAGS} )
  uniquely_add_flags( CMAKE_MODULE_LINKER_FLAGS ${LEDA_LINKER_FLAGS} )

 # Setup is done
  set ( LEDA_SETUP TRUE )

endif()

