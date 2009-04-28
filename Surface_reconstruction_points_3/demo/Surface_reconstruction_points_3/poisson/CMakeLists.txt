# This is the CMake script for compiling the Surface_reconstruction_points_3 MFC demo.

project( Poisson ) 

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true) 
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif(COMMAND cmake_policy)

# Require packages new or improved since CGAL 3.4
include_directories (BEFORE ../../../../Surface_mesher/include)
include_directories (BEFORE ../../../../Point_set_processing_3/include)

# Include this package's headers first
include_directories( BEFORE ./ ./include ../../../include )

# Find CGAL
find_package(CGAL)
if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
endif()

# This demo requires MFC 
if (MSVC)
  message( STATUS "MFC library found" )
  ADD_DEFINITIONS( "-D_AFXDLL" )
  SET(CMAKE_MFC_FLAG 2) 
endif()

# This demo requires OpenGL
find_package(OpenGL)
if(OPENGL_FOUND)
  message( STATUS "OpenGL include:     ${OPENGL_INCLUDE_DIR}" )
  include_directories ( ${OPENGL_INCLUDE_DIR} )

  message( STATUS "OpenGL definitions: ${OPENGL_DEFINITIONS}" )
  add_definitions( ${OPENGL_DEFINITIONS} )

  message( STATUS "OpenGL libraries:   ${OPENGL_LIBRARIES}" )
  link_libraries( ${OPENGL_LIBRARIES} )
endif()

# This demo requires TAUCS 
find_package(TAUCS)
if(TAUCS_FOUND)
  include( ${TAUCS_USE_FILE} )
endif()

if(CGAL_FOUND AND MSVC AND OPENGL_FOUND AND TAUCS_FOUND)

  # VisualC++ optimization for applications dealing with large data
  if (MSVC)
    # Use /EHa option to catch stack overflows.
    # Note: TAUCS needs a stack >= 2MB. CGAL default is 10MB.
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHa")

    # Allow Windows applications to use up to 3GB of RAM
    SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")

    # Turn off stupid VC++ warnings
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4311 /wd4800 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")

    # Speed up debugger
    SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_HAS_ITERATOR_DEBUGGING=0")

    # Print new compilation options
    message( STATUS "USING DEBUG CXXFLAGS   = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}'" )
    message( STATUS "USING DEBUG EXEFLAGS   = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}'" )
    message( STATUS "USING RELEASE CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}'" )
    message( STATUS "USING RELEASE EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}'" )
  endif()

  # Temporary debugging stuff
  #ADD_DEFINITIONS( "-DDEBUG_TRACE" )

  # Create Poisson executable
  ADD_EXECUTABLE(Poisson WIN32 ChildFrm.cpp DialogOptions.cpp director.cpp MainFrm.cpp Poisson.cpp PoissonDoc.cpp PoissonView.cpp Poisson.rc)
  
  # Link the executable to CGAL and third-party libraries
  if ( CGAL_AUTO_LINK_ENABLED )    
    target_link_libraries(Poisson ${CGAL_3RD_PARTY_LIBRARIES} )
  else()
    target_link_libraries(Poisson ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
  endif()

else (CGAL_FOUND AND MSVC AND OPENGL_FOUND AND TAUCS_FOUND)

  set(POISSON_MISSING_DEPS "")

  if(NOT CGAL_FOUND)
    set(POISSON_MISSING_DEPS "the CGAL library, ${POISSON_MISSING_DEPS}")
  endif()

  if(NOT MSVC)
    set(POISSON_MISSING_DEPS "MFC, ${POISSON_MISSING_DEPS}")
  endif()

  if(NOT OPENGL_FOUND)
    set(POISSON_MISSING_DEPS "OpenGL, ${POISSON_MISSING_DEPS}")
  endif()

  if(NOT TAUCS_FOUND)
    set(POISSON_MISSING_DEPS "TAUCS, ${POISSON_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This program requires ${POISSON_MISSING_DEPS}and will not be compiled.")

endif (CGAL_FOUND AND MSVC AND OPENGL_FOUND AND TAUCS_FOUND)
