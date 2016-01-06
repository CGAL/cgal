# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Arithmetic_kernel_test )

cmake_minimum_required(VERSION 2.8.11)



find_package(CGAL QUIET COMPONENTS Core GMP MPFR)

if ( CGAL_FOUND AND GMP_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_VersionUtils )
  
  get_dependency_version( GMP )

  IS_VERSION_LESS("${GMP_VERSION}" "4.2.0" _IS_GMP_VERSION_TO_LOW)
  
  include_directories(../../include)
  include_directories(include)

  include( CGAL_CreateSingleSourceCGALProgram )

  find_package( MPFI )

  if( MPFI_FOUND )
    # for the testsuite, the version of MPFI shipped with RS is used; this
    # version needs GMP>=4.2, so we require this dependency only here and
    # not in FindMPFI.cmake
    if( _IS_GMP_VERSION_TO_LOW )
      message( STATUS
        "MPFI tests need GMP>=4.2, some of the tests will not be compiled" )
    else( _IS_GMP_VERSION_TO_LOW )
      include( ${MPFI_USE_FILE} )
      create_single_source_cgal_program( "GMP_arithmetic_kernel.cpp" )      
    endif( _IS_GMP_VERSION_TO_LOW )
  else( MPFI_FOUND )
    message( STATUS
      "MPFI is not present, some of the tests will not be compiled." )
  endif( MPFI_FOUND )

  create_single_source_cgal_program( "Arithmetic_kernel.cpp" )
  create_single_source_cgal_program( "LEDA_arithmetic_kernel.cpp" )
  create_single_source_cgal_program( "CORE_arithmetic_kernel.cpp" )
  create_single_source_cgal_program( "Get_arithmetic_kernel.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

