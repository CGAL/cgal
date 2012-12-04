# This module setups the compiler for the MPFI library.
# It assumes that find_package(MPFI) was already called.

if( MPFI_FOUND AND NOT MPFI_SETUP )

  if (GMP_FOUND AND MPFR_FOUND) 

    message( STATUS "UseMPFI" )
    message( STATUS "MPFI include:      ${MPFI_INCLUDE_DIR}" )
    message( STATUS "MPFI definitions:  ${MPFI_DEFINITIONS}" )
    message( STATUS "MPFI libraries:    ${MPFI_LIBRARIES}" )

    try_run( MPFI_TEST_RESULT
             COMPILED_MPFI_TEST
             "${CMAKE_BINARY_DIR}"
             "${CGAL_MODULES_DIR}/test_MPFI.cpp"
             CMAKE_FLAGS
                "-DINCLUDE_DIRECTORIES:
                 STRING=${MPFI_INCLUDE_DIR};${GMP_INCLUDE_DIR};${MPFR_INCLUDE_DIR}"
               "-DLINK_LIBRARIES:
                 STRING=${MPFI_LIBRARIES};${GMP_LIBRARIES};${MPFR_LIBRARIES}"
                "-DLINK_DIRECTORIES:
                 STRING=${MPFI_LIBRARIES_DIR};${GMP_LIBRARIES_DIRS};${MPFR_LIBRARIES}"
             COMPILE_OUTPUT_VARIABLE MPFI_TEST_COMPILATION_OUTPUT
           )

    if( COMPILED_MPFI_TEST AND MPFI_TEST_RESULT EQUAL 0)
      include_directories( SYSTEM ${MPFI_INCLUDE_DIR} )
      link_directories( ${MPFI_LIBRARIES_DIR} )
      add_definitions( ${MPFI_DEFINITIONS} "-DCGAL_USE_MPFI" )
      link_libraries( ${MPFI_LIBRARIES} )
    else( COMPILED_MPFI_TEST AND MPFI_TEST_RESULT EQUAL 0)
      if (CGAL_ENABLE_PRECONFIG) 
        message( STATUS "MPFI is incorrectly configured with CGAL" )
      else()
        message( STATUS "MPFI is incorrectly configured on this system" )
      endif()
      message( STATUS
        "Output of the failed MPFI test was:\n${MPFI_TEST_COMPILATION_OUTPUT}" )
      message( STATUS "End of the MPFI test output" )
    endif( COMPILED_MPFI_TEST AND MPFI_TEST_RESULT EQUAL 0)

    set( MPFI_SETUP TRUE )

  else()

    message( STATUS "MPFI needs GMP and MPFR" )

  endif()

endif( MPFI_FOUND AND NOT MPFI_SETUP )

