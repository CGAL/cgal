# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Surface_mesh_parameterization_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

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

  include( ${CGAL_USE_FILE} )

  include( CreateSingleSourceCGALProgram )
  include( GeneratorSpecificSettings )

  # Link with BLAS, LAPACK and TAUCS (optional)
  include( CGAL_SetupTAUCS )

  # Link with Boost.ProgramOptions (optional)
  include( CGAL_SetupBoostProgramOptions )

  include_directories (BEFORE include)

  create_single_source_cgal_program( "Authalic_parameterization.cpp" )
  create_single_source_cgal_program( "Complete_parameterization_example.cpp" )
  create_single_source_cgal_program( "Mesh_cutting_parameterization.cpp" )
  create_single_source_cgal_program( "polyhedron_ex_parameterization.cpp" )
  create_single_source_cgal_program( "Simple_parameterization.cpp" )
  create_single_source_cgal_program( "Square_border_parameterization.cpp" )
  create_single_source_cgal_program( "Taucs_parameterization.cpp" )

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

