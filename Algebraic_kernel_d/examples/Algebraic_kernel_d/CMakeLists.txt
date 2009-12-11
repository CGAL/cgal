project( RS_examples )

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
    include_directories( BEFORE ../../include )
    create_single_source_cgal_program( "Construct_algebraic_real_1.cpp" )
    create_single_source_cgal_program( "Solve_1.cpp" )
    create_single_source_cgal_program( "Compare_1.cpp" )
    create_single_source_cgal_program( "Isolate_1.cpp" )
    create_single_source_cgal_program( "Sign_at_1.cpp" )
  else( RS_FOUND )
    message(STATUS "NOTICE: This program requires RS and will not be compiled.")
  endif( RS_FOUND )

else( CGAL_FOUND )

  message(STATUS
    "NOTICE: This program requires CGAL and GMP, and will not be compiled.")

endif( CGAL_FOUND )
