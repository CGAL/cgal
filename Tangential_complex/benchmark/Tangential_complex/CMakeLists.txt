# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Tangential_complex_benchmark )

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

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )
  
  find_package( TBB QUIET )

  if( TBB_FOUND )
    include(${TBB_USE_FILE})
    list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})
  endif()

  
  include( CGAL_CreateSingleSourceCGALProgram )

  find_package(Eigen3 3.1.0) 
  if (EIGEN3_FOUND) 
    include( ${EIGEN3_USE_FILE} ) 
    include_directories (BEFORE "../../include")
    include_directories (BEFORE "include")
    
    set (SOURCE_FILES "benchmark_tc.cpp")
    ADD_MSVC_PRECOMPILED_HEADER("StdAfx.h" "StdAfx.cpp" SOURCE_FILES)
    create_single_source_cgal_program( ${SOURCE_FILES} )
    
    set (SOURCE_FILES "benchmark_mesh_d.cpp")
    ADD_MSVC_PRECOMPILED_HEADER("StdAfx.h" "StdAfx.cpp" SOURCE_FILES)
    create_single_source_cgal_program( ${SOURCE_FILES} )

  else()
    message(STATUS "NOTICE: Some of the executables in this directory need Eigen 3.1 (or greater) and will not be compiled.")
  endif()

else()
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()

