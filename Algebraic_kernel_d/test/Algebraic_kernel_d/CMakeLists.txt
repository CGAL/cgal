project( RS_test )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
endif()

find_package( CGAL QUIET )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
  find_package( RS )
  if( RS_FOUND )
    include( ${RS_USE_FILE} )
    include_directories (BEFORE ../../include)
    create_single_source_cgal_program( "generic_test_integer.cpp" )
    create_single_source_cgal_program( "generic_test_rational.cpp" )
    create_single_source_cgal_program( "io_test.cpp" )
    create_single_source_cgal_program( "solve_integer.cpp" )
    create_single_source_cgal_program( "solve_rational.cpp" )
  else( RS_FOUND )
    message(STATUS "NOTICE: This program requires RS and will not be compiled.")
  endif( RS_FOUND )

else( CGAL_FOUND )

  message(STATUS
    "NOTICE: This program requires CGAL and GMP, and will not be compiled.")

endif( CGAL_FOUND )
