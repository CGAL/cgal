project( Scale_space )

# Check CMake version
cmake_minimum_required( VERSION 2.6.2 )
if( "${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6 )
  if( "${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3 )
    cmake_policy( VERSION 2.8.4 )
  else()
    cmake_policy( VERSION 2.6 )
  endif()
endif()


# Preprocessor directives
if( WIN32 )
  add_definitions( "/W3 /D_CRT_SECURE_NO_WARNINGS /wd4005 /wd4307 /wd4308 /nologo" )
endif()
if( UNIX )
  add_definitions( "-fpermissive" )
endif()

# Enable openMP for parallel processing
if( WIN32 )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp" )
  set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /openmp" )
endif()


find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include" $ENV{EIGEN3_DIR}/ )

  create_single_source_cgal_program( "scale_space.cpp" )

else()
  
  message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()
