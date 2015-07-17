# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Mesh_3_benchmark )

cmake_minimum_required(VERSION 2.8.11)


# Creates a new CMake option, turned ON by default
option(ACTIVATE_MSVC_PRECOMPILED_HEADERS
  "Activate precompiled headers in MSVC"
  OFF)

# Macro to add precompiled headers for MSVC
# This function does two things:
# 1. Enable precompiled headers on each file which is listed in "SourcesVar".
# 2. Add the content of "PrecompiledSource" (e.g. "StdAfx.cpp") to "SourcesVar".
MACRO(ADD_MSVC_PRECOMPILED_HEADER PrecompiledHeader PrecompiledSource SourcesVar)
  IF(MSVC AND ACTIVATE_MSVC_PRECOMPILED_HEADERS)
    GET_FILENAME_COMPONENT(PrecompiledBasename ${PrecompiledHeader} NAME_WE)
    SET(Sources ${${SourcesVar}})

    SET_SOURCE_FILES_PROPERTIES(${PrecompiledSource}
                                PROPERTIES COMPILE_FLAGS "/Yc\"${PrecompiledHeader}\"")
    SET_SOURCE_FILES_PROPERTIES(${Sources}
                                PROPERTIES COMPILE_FLAGS "/Yu\"${PrecompiledHeaders}\" /FI\"${PrecompiledHeader}\"")
    # Add precompiled header to SourcesVar
    LIST(APPEND ${SourcesVar} ${PrecompiledSource})
  ENDIF(MSVC AND ACTIVATE_MSVC_PRECOMPILED_HEADERS)
ENDMACRO(ADD_MSVC_PRECOMPILED_HEADER)
# The compiler might need more memory because of precompiled headers
if(MSVC AND ACTIVATE_MSVC_PRECOMPILED_HEADERS AND NOT(MSVC_VERSION LESS 1310))
  set(CGAL_C_FLAGS "${CGAL_C_FLAGS} /Zm1000")
  set(CGAL_CXX_FLAGS "${CGAL_CXX_FLAGS} /Zm1000")
endif()

include_directories(../../include)
include_directories(../../../Triangulation_3/include)
include_directories(../../../STL_Extension/include)
include_directories(../../../AABB_tree/include)
add_definitions(-DCGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
                -DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)

if ( MESH_3_VERBOSE )
  add_definitions(-DCGAL_MESH_3_VERBOSE)
endif()

find_package(CGAL COMPONENTS ImageIO)

if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )

  # Activate concurrency ? (turned OFF by default)
  option(ACTIVATE_CONCURRENT_MESH_3
    "Activate parallelism in Mesh_3"
    OFF)
    
  # And add -DCGAL_CONCURRENT_MESH_3 if that option is ON
  if( ACTIVATE_CONCURRENT_MESH_3 )
    add_definitions( -DCGAL_CONCURRENT_MESH_3 )
    find_package( TBB REQUIRED )	
  else()
    option( LINK_WITH_TBB
      "Link with TBB anyway so we can use TBB timers for profiling"
      ON)
    if( LINK_WITH_TBB )
      find_package( TBB )
    endif( LINK_WITH_TBB )
  endif()
  
  if( TBB_FOUND )
    include(${TBB_USE_FILE})
    list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})
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
    list(APPEND CGAL_3RD_PARTY_LIBRARIES ${Boost_LIBRARIES})
  endif()

  if ( Boost_FOUND AND Boost_VERSION GREATER 103400 )
    include( CGAL_CreateSingleSourceCGALProgram )

    # Compilable benchmark
    set (BENCHMARK_SOURCE_FILES "concurrency.cpp")
    ADD_MSVC_PRECOMPILED_HEADER("StdAfx.h" "StdAfx.cpp" BENCHMARK_SOURCE_FILES)
    create_single_source_cgal_program( ${BENCHMARK_SOURCE_FILES} )
    
  else()
    message(STATUS "NOTICE: This program requires Boost >= 1.34.1, and will not be compiled.")
  endif()

else()
  message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()

