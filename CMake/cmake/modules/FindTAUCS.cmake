# Try to find the TAUCS libraries
# TAUCS_FOUND - system has TAUCS lib
# TAUCS_INCLUDE_DIR - the TAUCS include directory
# TAUCS_LIBRARIES_DIR - Directory where the TAUCS libraries are located
# TAUCS_LIBRARIES - the TAUCS libraries
# TAUCS_IN_AUXILIARY - TRUE if the TAUCS found is the one distributed with CGAL in the auxiliary folder

# TODO: support MacOSX

include(FindPackageHandleStandardArgs)
include(GeneratorSpecificSettings)

# Is it already configured?
if (TAUCS_INCLUDE_DIR AND TAUCS_LIBRARIES ) 
   
  set(TAUCS_FOUND TRUE)
  
else()  

  # Look first for the TAUCS distributed with CGAL in auxiliary/taucs
  find_path(TAUCS_INCLUDE_DIR 
            NAMES taucs.h 
            PATHS ${CGAL_SOURCE_DIR}/auxiliary/taucs/include
                  ENV TAUCS_INC_DIR
  	        DOC "The directory containing the TAUCS header files"
           )

  if ( TAUCS_INCLUDE_DIR ) 

     if ( TAUCS_INCLUDE_DIR STREQUAL "${CGAL_SOURCE_DIR}/auxiliary/taucs/include" )
       set( TAUCS_IN_CGAL_AUXILIARY TRUE )
       set( TAUCS_LIB_SEARCH_PATHS ${CGAL_SOURCE_DIR}/auxiliary/taucs/lib )
     endif()
     
    if ( AUTO_LINK_ENABLED )
      set( TAUCS_NAMES "libtaucs${TOOLSET}-mt" "libtaucs${TOOLSET}-mt-gd" "libtaucs${TOOLSET}-mt-o" "libtaucs${TOOLSET}-mt-g" )
    else()
      set( TAUCS_NAMES "taucs" )
    endif()  
    
    find_library(TAUCS_LIBRARIES 
                 NAMES ${TAUCS_NAMES} 
                 PATHS ${TAUCS_LIB_SEARCH_PATHS}
                       ENV TAUCS_LIB_DIR
                 DOC "Path to the TAUCS library"
                )
     
    if ( TAUCS_LIBRARIES ) 
      get_filename_component(TAUCS_LIBRARIES_DIR ${TAUCS_LIBRARIES} PATH)
    endif()
    
  endif()  
  
  # Attempt to load a user-defined configuration for TAUCS if couldn't be found
  if ( NOT TAUCS_INCLUDE_DIR OR NOT TAUCS_LIBRARIES )
    include( TAUCSConfig OPTIONAL )
  endif()
  
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(TAUCS "DEFAULT_MSG" TAUCS_INCLUDE_DIR TAUCS_LIBRARIES)
  
endif()


