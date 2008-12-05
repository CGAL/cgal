# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Surface_mesh_parameterization_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
endif()

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  # Optimization for applications dealing with large data
  if (MSVC)
    # Use /EHa option to catch stack overflows.
    # Note: TAUCS needs a stack >= 2MB. CGAL default is 10MB.
    SET (CGAL_CXX_FLAGS_INIT "/DWIN32 /D_WINDOWS /W2 /Zm1000 /EHa /GR")

    # Allow Windows applications to use up to 3GB of RAM
    SET (CGAL_EXE_LINKER_FLAGS_INIT "${CGAL_EXE_LINKER_FLAGS_INIT} /LARGEADDRESSAWARE")

    # Turn off stupid VC++ warnings
    SET (CGAL_CXX_FLAGS_INIT "${CGAL_CXX_FLAGS_INIT} /wd4267 /wd4311 /wd4800 /wd4503 /wd4244 /wd4345")

    # Speed up VC++ debugger
    SET (CGAL_CXX_FLAGS_DEBUG_INIT "${CGAL_CXX_FLAGS_DEBUG_INIT} /D_HAS_ITERATOR_DEBUGGING=0")
  endif()

  # g++ 4.1.2 optimizer creates an infinite loop in Surface_mesh_parameterization if -O2 or -O3
  if (CMAKE_COMPILER_IS_GNUCXX)
    SET (CGAL_CXX_FLAGS_RELEASE_INIT "${CGAL_CXX_FLAGS_RELEASE_INIT} -O")
  endif()

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  # Link with Boost.ProgramOptions (optional)
  find_package(Boost QUIET COMPONENTS program_options)
  if(Boost_PROGRAM_OPTIONS_FOUND)
    if( CGAL_AUTO_LINK_ENABLED )
      message( STATUS "Boost.ProgramOptions library: found" )
    else()
      message( STATUS "Boost.ProgramOptions library: ${Boost_PROGRAM_OPTIONS_LIBRARY}" )
    endif()
    add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
    if ( NOT CGAL_AUTO_LINK_ENABLED )
      link_libraries( ${Boost_PROGRAM_OPTIONS_LIBRARY} )
    endif()
  endif()

  include_directories (BEFORE include)

  # examples that do not use TAUCS
  create_single_source_cgal_program( "Authalic_parameterization.cpp" )
  create_single_source_cgal_program( "Mesh_cutting_parameterization.cpp" )
  create_single_source_cgal_program( "Simple_parameterization.cpp" )
  create_single_source_cgal_program( "Square_border_parameterization.cpp" )

  # Link with BLAS, LAPACK and TAUCS (optional)
  find_package(TAUCS)
  if(TAUCS_FOUND)

    include( ${TAUCS_USE_FILE} )

    # examples that use TAUCS
    create_single_source_cgal_program( "Complete_parameterization_example.cpp" )
    create_single_source_cgal_program( "Taucs_parameterization.cpp" )

  else(TAUCS_FOUND)

    message(STATUS "NOTICE: This directory requires the TAUCS library, and some examples will not be compiled.")

  endif(TAUCS_FOUND)

  # example that can use TAUCS or not, depending on the macro CGAL_USE_TAUCS
  create_single_source_cgal_program( "polyhedron_ex_parameterization.cpp" )

else()

  message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

