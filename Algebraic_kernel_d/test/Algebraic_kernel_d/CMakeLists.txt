# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Algebraic_kernel_d_test )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
endif()

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
  include( CGAL_VersionUtils )
  find_package( MPFI )

  if( MPFI_FOUND )
    include( ${MPFI_USE_FILE} )
  endif( MPFI_FOUND )

  find_package( RS )

  if( RS_FOUND )
    include( ${RS_USE_FILE} )
  endif( RS_FOUND )

  include_directories (BEFORE ../../include)
  include_directories (BEFORE include)

  create_single_source_cgal_program( "Algebraic_curve_kernel_2.cpp" )
  create_single_source_cgal_program( "algebraic_curve_kernel_2_tools.cpp" )
  create_single_source_cgal_program( "Algebraic_kernel_d_1_LEDA.cpp" )
  create_single_source_cgal_program( "Algebraic_kernel_d_1_CORE.cpp" )
  create_single_source_cgal_program( "Algebraic_kernel_d_1_GMP.cpp" )
  create_single_source_cgal_program( "Algebraic_kernel_d_2.cpp" )
  create_single_source_cgal_program( "Algebraic_kernel_rs_gmpq_d_1.cpp" )
  create_single_source_cgal_program( "Algebraic_kernel_rs_gmpz_d_1.cpp" )
  create_single_source_cgal_program( "Algebraic_real_d_1.cpp" )
  create_single_source_cgal_program( "Bitstream_descartes.cpp" )
  create_single_source_cgal_program( "Curve_analysis_2.cpp" )
  create_single_source_cgal_program( "Curve_pair_analysis_2.cpp" )
  create_single_source_cgal_program( "Descartes.cpp" )
  create_single_source_cgal_program( "Real_embeddable_traits_extension.cpp" )

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

