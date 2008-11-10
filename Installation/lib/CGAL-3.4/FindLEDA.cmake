if (LEDA_INCLUDE_DIR AND LEDA_LIBRARIES ) 
   
  set(LEDA_FOUND TRUE)
  
else()  
  
  find_path(LEDA_INCLUDE_DIR 
            NAMES "LEDA/basic.h" "LEDA/system/basic.h"
            PATHS ENV LEDA_INC_DIR
            DOC "The directory containing the LEDA header files WITHOUT the LEDA prefix"
          )
          
  find_library(LEDA_LIBRARY_RELEASE NAMES "leda"
                PATHS ENV LEDA_LIB_DIR
                DOC "Path to the LEDA library"
              )
              
  find_library(LEDA_LIBRARY_DEBUG NAMES "ledaD"
                PATHS ENV LEDA_LIB_DIR
                DOC "Path to the LEDA library"
              )
              
  if ( NOT LEDA_INCLUDE_DIR )
    set( LEDA_INCLUDE_DIR   "$ENV{LEDA_INC_DIR}" CACHE FILEPATH "The directory containing the LEDA header files WITHOUT the LEDA prefix" FORCE )
  endif()
    
  if ( NOT LEDA_LIBRARY_RELEASE )
    set( LEDA_LIBRARY_RELEASE  "$ENV{LEDA_LIBRARY_RELEASE}"  CACHE FILEPATH "The LEDA release-mode libraries" FORCE )
  endif()

  if ( NOT LEDA_LIBRARY_DEBUG )
    set( LEDA_LIBRARY_DEBUG "$ENV{LEDA_LIBRARY_DEBUG}"  CACHE FILEPATH "The LEDA debug-mode libraries" FORCE )
  endif()
  
  if ( "${CMAKE_BUILD_TYPE}" STREQUAL "Release" )
    if ( LEDA_LIBRARY_RELEASE )
      set( LEDA_LIBRARIES "${LEDA_LIBRARY_RELEASE}" )
    endif()  
  else()
    if ( LEDA_LIBRARY_DEBUG )
      set( LEDA_LIBRARIES "${LEDA_LIBRARY_DEBUG}" )
    endif()  
  endif()
  
  if ( NOT LEDA_DEFINITIONS )
    set( LEDA_DEFINITIONS "$ENV{LEDA_DEFINITIONS}" CACHE STRING   "Definitions for the LEDA library" FORCE )
  endif()  
  
  if ( NOT LEDA_CXX_FLAGS )
    set( LEDA_CXX_FLAGS "$ENV{LEDA_CXX_FLAGS}" CACHE STRING   "Compiler flags for the LEDA library" FORCE )
  endif()  
  
  if ( NOT LEDA_LINKER_FLAGS )
    set( LEDA_LINKER_FLAGS "$ENV{LEDA_LINKER_FLAGS}" CACHE STRING   "Linker flags for the LEDA library" FORCE )
  endif()  
  
  if ( NOT "$ENV{LEDA_VERSION}" STREQUAL "" )
    set( LEDA_VERSION "$ENV{LEDA_VERSION}" CACHE STRING "The version of the LEDA library" FORCE )
  endif()
  
  if ( LEDA_VERSION )
    if ( NOT "${LEDA_DEFINITIONS}" MATCHES "-DCGAL_LEDA_VERSION=${LEDA_VERSION}" )
      set( LEDA_DEFINITIONS "${LEDA_DEFINITIONS}" "-DCGAL_LEDA_VERSION=${LEDA_VERSION}" CACHE STRING "Definitions for the LEDA library" FORCE )
    endif()
  endif()
  
  if ( LEDA_INCLUDE_DIR AND LEDA_LIBRARIES)
    set(LEDA_FOUND TRUE)
  endif()
  
endif()
