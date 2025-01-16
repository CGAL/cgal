# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Arithmetic_kernel_Tests)

find_package(CGAL REQUIRED COMPONENTS Core)

find_package(GMP)

if(GMP_FOUND)

  include(CGAL_VersionUtils)

  get_dependency_version(GMP)

  is_version_less("${GMP_VERSION}" "4.2.0" _IS_GMP_VERSION_TO_LOW)

  include_directories(include)

  find_package(MPFI QUIET)

  if(MPFI_FOUND)
    include(${MPFI_USE_FILE})

    # for the testsuite, the version of MPFI shipped with RS is used; this
    # version needs GMP>=4.2, so we require this dependency only here and
    # not in FindMPFI.cmake
    if(_IS_GMP_VERSION_TO_LOW)
      message(STATUS "NOTICE: MPFI tests need GMP>=4.2, some of the tests will not be compiled")
    else(_IS_GMP_VERSION_TO_LOW)
      include(${MPFI_USE_FILE})
      create_single_source_cgal_program("GMP_arithmetic_kernel.cpp")
    endif(_IS_GMP_VERSION_TO_LOW)
  else(MPFI_FOUND)
    message(STATUS "NOTICE: MPFI is not present, some of the tests will not be compiled.")
  endif(MPFI_FOUND)

  create_single_source_cgal_program("Arithmetic_kernel.cpp")
  create_single_source_cgal_program("LEDA_arithmetic_kernel.cpp")
  create_single_source_cgal_program("CORE_arithmetic_kernel.cpp")
  create_single_source_cgal_program("GMPXX_arithmetic_kernel.cpp")
  create_single_source_cgal_program("Get_arithmetic_kernel.cpp")

else()

  message("NOTICE: This project requires GMP support, and will not be compiled.")

endif()
