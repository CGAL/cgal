# This module setups the compiler for the MPFI library.
# It assumes that find_package(MPFI) was already called.

if( NOT CGAL_MPFI_SETUP )

  if( MPFI_FOUND )

    message( STATUS "MPFI include:      ${MPFI_INCLUDE_DIR}" )
    message( STATUS "MPFI libraries:    ${MPFI_LIBRARIES}" )
    message( STATUS "MPFI definitions:  ${MPFI_DEFINITIONS}" )

    try_run( MPFI_TEST_RESULT
             COMPILED_MPFI_TEST
             "${CMAKE_BINARY_DIR}"
             "${CGAL_MODULES_DIR}/test_MPFI.cpp"
             CMAKE_FLAGS
               "-DINCLUDE_DIRECTORIES:
                 STRING=${MPFI_INCLUDE_DIR};${CGAL_3RD_PARTY_INCLUDE_DIRS}"
               "-DLINK_LIBRARIES:
                 STRING=${MPFI_LIBRARIES};${CGAL_3RD_PARTY_LIBRARIES}"
               "-DLINK_DIRECTORIES:
                 STRING=${MPFI_LIBRARIES_DIR};${CGAL_3RD_PARTY_LIBRARIES_DIRS}"
             COMPILE_OUTPUT_VARIABLE MPFI_TEST_COMPILATION_OUTPUT
           )

    if( COMPILED_MPFI_TEST AND MPFI_TEST_RESULT EQUAL 0)
      include_directories( ${MPFI_INCLUDE_DIR} )
      link_directories( ${MPFI_LIBRARIES_DIR} )
      add_definitions( ${MPFI_DEFINITIONS} "-DCGAL_USE_MPFI" )
      link_libraries( ${MPFI_LIBRARIES} )
    else( COMPILED_MPFI_TEST AND MPFI_TEST_RESULT EQUAL 0)
      message( STATUS "MPFI was incorrectly configured on this system" )
      message( STATUS
        "Output of the failed MPFI test was:\n${MPFI_TEST_COMPILATION_OUTPUT}" )
      message( STATUS "End of the MPFI test output" )
    endif( COMPILED_MPFI_TEST AND MPFI_TEST_RESULT EQUAL 0)

  endif( MPFI_FOUND )

  set( CGAL_MPFI_SETUP TRUE )

endif( NOT CGAL_MPFI_SETUP )
