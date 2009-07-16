# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Surface_mesh_parameterization_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
endif()
 
# Find CGAL
find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )

  # VisualC++ optimization for applications dealing with large data
  if (MSVC)
    # Use /EHa option to catch stack overflows.
    # Note: TAUCS needs a stack >= 2MB. CGAL default is 10MB.
    string(REGEX REPLACE "/EH[asc]*" "/EHa" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

    # Allow Windows applications to use up to 3GB of RAM
    SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")

    # Turn off stupid VC++ warnings
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4311 /wd4800 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")
    
    # Print new compilation options
    message( STATUS "USING DEBUG CXXFLAGS   = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}'" )
    message( STATUS "USING DEBUG EXEFLAGS   = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}'" )
    message( STATUS "USING RELEASE CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}'" )
    message( STATUS "USING RELEASE EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}'" )
  endif()

  # g++ 4.1 creates an infinite loop in Surface_mesh_parameterization
  # if strict-aliasing optimization is on (implied by -O2, -O3 and -Os).
  if (CMAKE_COMPILER_IS_GNUCXX)
    SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-strict-aliasing")
    SET (CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   -fno-strict-aliasing")

    # Print new compilation options
    string( TOUPPER "${CMAKE_BUILD_TYPE}" CGAL_BUILD_TYPE_UPPER )
    message( STATUS "USING CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CGAL_BUILD_TYPE_UPPER}}'" )
    message( STATUS "USING EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_${CGAL_BUILD_TYPE_UPPER}}'" )
  endif()

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

  include_directories (BEFORE ../../include)
  include_directories (BEFORE include)

  # Executables that do *not* use TAUCS
  create_single_source_cgal_program( "Authalic_parameterization.cpp" )
  create_single_source_cgal_program( "Mesh_cutting_parameterization.cpp" )
  create_single_source_cgal_program( "Simple_parameterization.cpp" )
  create_single_source_cgal_program( "Square_border_parameterization.cpp" )

  # Link with BLAS, LAPACK and TAUCS (optional)
  find_package(TAUCS)
  if(TAUCS_FOUND)

    include( ${TAUCS_USE_FILE} )

    # Executables that require TAUCS (thus BLAS and LAPACK)
    create_single_source_cgal_program( "Complete_parameterization_example.cpp" )
    create_single_source_cgal_program( "Taucs_parameterization.cpp" )

  else(TAUCS_FOUND)

    message(STATUS "NOTICE: This directory requires the TAUCS library, and some features will not be available.")

  endif(TAUCS_FOUND)

  # Executable that can use TAUCS or not, depending on the macro CGAL_USE_TAUCS
  create_single_source_cgal_program( "polyhedron_ex_parameterization.cpp" )

else()

  message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

