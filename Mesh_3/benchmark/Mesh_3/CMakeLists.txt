# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Mesh_3_benchmark )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

# Creates a new CMake option, turned ON by default
option(ACTIVATE_MSVC_PRECOMPILED_HEADERS
  "Activate precompiled headers in MSVC"
  ON)
  
# Macro to add precompiled headers for MSVC
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

include_directories(../../include)
include_directories(../../../Triangulation_3/include)
include_directories(../../../STL_Extension/include)
include_directories(../../../AABB_tree/include)
add_definitions(-DCGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
                -DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)

# Disable warning 4503 in MSVC ("The decorated name was longer than the 
# compiler limit and was truncated")
if(MSVC)
  add_definitions( "/wd4503" )
 endif()

# Add specific Find.cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../../../Installation/cmake/modules/ )

if ( MESH_3_VERBOSE )
  add_definitions(-DCGAL_MESH_3_VERBOSE)
endif()

find_package(CGAL COMPONENTS ImageIO)

if ( CGAL_FOUND )

  # Activate concurrency ? (turned OFF by default)
  option(ACTIVATE_CONCURRENT_MESH_3
    "Activate parallelism in Mesh_3"
    OFF)
    
  # And add -DCONCURRENT_MESH_3 if that option is ON
  if( ACTIVATE_CONCURRENT_MESH_3 )
    add_definitions( -DCONCURRENT_MESH_3 )
    find_package( TBB REQUIRED )
  else( ACTIVATE_CONCURRENT_MESH_3 )
    option( LINK_WITH_TBB
      "Link with TBB anyway so we can use TBB timers for profiling"
      ON)
    if( LINK_WITH_TBB )
      find_package( TBB )
    endif( LINK_WITH_TBB )
  endif()
  
  if( TBB_FOUND )
    add_definitions( -DNOMINMAX -DLINKED_WITH_TBB )
    include_directories ( ${TBB_INCLUDE_DIRS} )
    link_directories( ${TBB_LIBRARY_DIRS} )
  endif( TBB_FOUND )
  
  include( ${CGAL_USE_FILE} )
  find_package(Boost)

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

