# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( Algebraic_kernel_d_test )

cmake_minimum_required(VERSION 2.8.11)


# CGAL and its components
find_package( CGAL QUIET COMPONENTS Core RS3)

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )

include( CGAL_CreateSingleSourceCGALProgram )

# Boost and its components
find_package( Boost )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory
include_directories( BEFORE include )

# include for local package
include_directories( BEFORE ../../include )

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "cyclic.cpp" )
create_single_source_cgal_program( "Algebraic_curve_kernel_2.cpp" )
create_single_source_cgal_program( "algebraic_curve_kernel_2_tools.cpp" )
create_single_source_cgal_program( "Algebraic_kernel_d_1_LEDA.cpp" )

create_single_source_cgal_program( "Algebraic_kernel_d_1_CORE_Integer_rational.cpp" )
create_single_source_cgal_program( "Algebraic_kernel_d_1_CORE_SqrtII_rational.cpp" )
create_single_source_cgal_program( "Algebraic_kernel_d_1_CORE_SqrtRI_rational.cpp" )
create_single_source_cgal_program( "Algebraic_kernel_d_1_CORE_SqrtRR_rational.cpp" )

create_single_source_cgal_program( "Algebraic_kernel_d_1_GMP.cpp" )
create_single_source_cgal_program( "Algebraic_kernel_d_2.cpp" )
create_single_source_cgal_program( "Algebraic_real_d_1.cpp" )
create_single_source_cgal_program( "Bitstream_descartes.cpp" )
create_single_source_cgal_program( "Curve_analysis_2.cpp" )
create_single_source_cgal_program( "Curve_pair_analysis_2.cpp" )
create_single_source_cgal_program( "Descartes.cpp" )
create_single_source_cgal_program( "Real_embeddable_traits_extension.cpp" )
if(RS_FOUND)
  create_single_source_cgal_program( "Algebraic_kernel_rs_gmpq_d_1.cpp" )
  create_single_source_cgal_program( "Algebraic_kernel_rs_gmpz_d_1.cpp" )
else()
  message(STATUS "NOTICE: Some tests require the RS library, and will not be compiled.")
endif()
