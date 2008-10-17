if (LEDA_INCLUDE_DIR AND LEDA_LIBRARIES ) 
   
  set(LEDA_FOUND TRUE)
  
else()  
  
  find_path(LEDA_INCLUDE_DIR 
            NAMES "LEDA/basic.h" "LEDA/system/basic.h"
            PATHS ENV LEDA_INC_DIR
            DOC "The directory containing the LEDA header files WITHOUT the LEDA prefix"
          )
          
  find_library(LEDA_LIBRARIES NAMES "leda"
                PATHS ENV LEDA_LIB_DIR
                DOC "Path to the LEDA library"
              )
              
  if ( NOT LEDA_INCLUDE_DIR )
    set( LEDA_INCLUDE_DIR   "$ENV{LEDA_INC_DIR}" CACHE FILEPATH "The directory containing the LEDA header files WITHOUT the LEDA prefix" )
  endif()
    
  if ( NOT LEDA_LIBRARIES )
    set( LEDA_LIBRARIES  "$ENV{LEDA_LIBRARIES}"  CACHE FILEPATH "The LEDA libraries" )
  endif()

  if ( NOT LEDA_DEFINITIONS )
    set( LEDA_DEFINITIONS "$ENV{LEDA_DEFINITIONS}" CACHE STRING   "Definitions for the LEDA library" )
  endif()  
  
  if ( LEDA_INCLUDE_DIR AND LEDA_LIBRARIES)
    set(LEDA_FOUND TRUE)
  endif()
  
endif()
