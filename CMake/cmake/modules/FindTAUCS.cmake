# Try to find the TAUCS libraries
# TAUCS_FOUND - system has TAUCS lib
# TAUCS_INCLUDE_DIR - the TAUCS include directory
# TAUCS_LIBRARIES_DIR - Directory where the TAUCS libraries are located
# TAUCS_LIBRARIES - the TAUCS libraries
# TAUCS_IN_CGAL_AUXILIARY - TRUE if the TAUCS found is the one distributed with CGAL in the auxiliary folder

# TODO: support MacOSX

include(FindPackageHandleStandardArgs)
include(GeneratorSpecificSettings)

# Is it already configured?
if (TAUCS_INCLUDE_DIR AND TAUCS_LIBRARIES_DIR ) 
   
  set(TAUCS_FOUND TRUE)
  
else()  

  # Look first for the TAUCS distributed with CGAL in auxiliary/taucs
  find_path(TAUCS_INCLUDE_DIR 
            NAMES taucs.h 
            PATHS ${CGAL_SOURCE_DIR}/auxiliary/taucs/include
                  ENV TAUCS_INC_DIR
  	        DOC "The directory containing the TAUCS header files"
           )

  if ( TAUCS_INCLUDE_DIR STREQUAL "${CGAL_SOURCE_DIR}/auxiliary/taucs/include" )
    set( TAUCS_IN_CGAL_AUXILIARY TRUE CACHE INTERNAL "" )
  endif()
  
  if ( AUTO_LINK_ENABLED )
  
    find_path(TAUCS_LIBRARIES_DIR 
              NAMES "libtaucs${TOOLSET}-mt" "libtaucs${TOOLSET}-mt-gd" "libtaucs${TOOLSET}-mt-o" "libtaucs${TOOLSET}-mt-g"
              PATHS ${CGAL_SOURCE_DIR}/auxiliary/taucs/lib
                    ENV TAUCS_LIB_DIR
              DOC "Directory containing the TAUCS library"
             ) 
    
  else()
  
    find_library(TAUCS_LIBRARIES NAMES "taucs"
                 PATHS ENV TAUCS_LIB_DIR
                 DOC "Path to the TAUCS library"
                )
                
    if ( TAUCS_LIBRARIES ) 
      get_filename_component(TAUCS_LIBRARIES_DIR ${TAUCS_LIBRARIES} PATH CACHE)
    endif()
    
  endif()  
  
  # Attempt to load a user-defined configuration for TAUCS if couldn't be found
  if ( NOT TAUCS_INCLUDE_DIR OR NOT TAUCS_LIBRARIES_DIR )
    include( TAUCSConfig OPTIONAL )
  endif()
  
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(TAUCS "DEFAULT_MSG" TAUCS_INCLUDE_DIR TAUCS_LIBRARIES_DIR)
  
endif()

mark_as_advanced(TAUCS_INCLUDE_DIR)
mark_as_advanced(TAUCS_LIBRARIES)
mark_as_advanced(TAUCS_LIBRARIES_DIR)
