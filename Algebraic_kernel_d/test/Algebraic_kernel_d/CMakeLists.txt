cmake_minimum_required(VERSION 3.1...3.23)
project(Algebraic_kernel_d_Tests)

# CGAL and its components
find_package(CGAL REQUIRED COMPONENTS Core)

find_package(MPFI QUIET)
if(MPFI_FOUND)
  message(STATUS "Found MPFI")
  include(${MPFI_USE_FILE})
endif()

find_package(RS3 QUIET)
if(RS3_FOUND)
  message(STATUS "Found RS3")
  include(${RS3_USE_FILE})
endif()

# include for local directory
include_directories(BEFORE include)

create_single_source_cgal_program("cyclic.cpp")
create_single_source_cgal_program("Descartes.cpp")

if(NOT CGAL_DISABLE_GMP)
  create_single_source_cgal_program("Algebraic_curve_kernel_2.cpp")
  create_single_source_cgal_program("algebraic_curve_kernel_2_tools.cpp")
  create_single_source_cgal_program("Algebraic_kernel_d_1_LEDA.cpp")

  create_single_source_cgal_program("Algebraic_kernel_d_1_CORE_Integer_rational.cpp")
  create_single_source_cgal_program("Algebraic_kernel_d_1_CORE_SqrtII_rational.cpp")
  create_single_source_cgal_program("Algebraic_kernel_d_1_CORE_SqrtRI_rational.cpp")
  create_single_source_cgal_program("Algebraic_kernel_d_1_CORE_SqrtRR_rational.cpp")

  create_single_source_cgal_program("Algebraic_kernel_d_1_GMP.cpp")
  create_single_source_cgal_program("Algebraic_kernel_d_2.cpp")
  create_single_source_cgal_program("Algebraic_real_d_1.cpp")
  create_single_source_cgal_program("Bitstream_descartes.cpp")
  create_single_source_cgal_program("Curve_analysis_2.cpp")
  create_single_source_cgal_program("Curve_pair_analysis_2.cpp")
  create_single_source_cgal_program("Real_embeddable_traits_extension.cpp")

  if(RS_FOUND)
    create_single_source_cgal_program("Algebraic_kernel_rs_gmpq_d_1.cpp")
    create_single_source_cgal_program("Algebraic_kernel_rs_gmpz_d_1.cpp")
  else()
    message(STATUS "NOTICE: Some tests require the RS library, and will not be compiled.")
  endif()
else()
  message(STATUS "NOTICE: Some tests require GMP support, and will not be compiled.")
endif()
