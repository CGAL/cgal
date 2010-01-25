project( RS_examples )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
endif()

find_package( CGAL QUIET )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_VersionUtils )
  include( CGAL_CreateSingleSourceCGALProgram )
  find_package( RS )
  if( RS_FOUND )
    IS_VERSION_LESS("${CGAL_GMP_VERSION}" "4.2.0" _IS_GMP_VERSION_TO_LOW)
    if(_IS_GMP_VERSION_TO_LOW)
      message(STATUS "NOTICE: These program uses the CGAL Algebraic Kernel, and its implementation using RS needs GMP>=4.2. Your GMP version is ${CGAL_GMP_VERSION}.")
    else()
      include( ${RS_USE_FILE} )
      include_directories( BEFORE ../../include )
      create_single_source_cgal_program( "Construct_algebraic_real_1.cpp" )
      create_single_source_cgal_program( "Solve_1.cpp" )
      create_single_source_cgal_program( "Compare_1.cpp" )
      create_single_source_cgal_program( "Isolate_1.cpp" )
      create_single_source_cgal_program( "Sign_at_1.cpp" )
    endif()
  else( RS_FOUND )
    message(STATUS "NOTICE: This program requires RS and will not be compiled.")
  endif( RS_FOUND )

else( CGAL_FOUND )

  message(STATUS
    "NOTICE: This program requires CGAL and GMP, and will not be compiled.")

endif( CGAL_FOUND )
