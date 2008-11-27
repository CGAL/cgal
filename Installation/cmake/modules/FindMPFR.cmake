# Try to find the MPFR libraries
# MPFR_FOUND - system has MPFR lib
# MPFR_INCLUDE_DIR - the MPFR include directory
# MPFR_LIBRARIES_DIR - Directory where the MPFR libraries are located
# MPFR_LIBRARIES - the MPFR libraries
# MPFR_IN_CGAL_AUXILIARY - TRUE if the MPFR found is the one distributed with CGAL in the auxiliary folder

# TODO: support MacOSX

include(CGAL_FindPackageHandleStandardArgs)
include(CGAL_GeneratorSpecificSettings)

# Is it already configured?
if (MPFR_INCLUDE_DIR AND MPFR_LIBRARIES_DIR ) 
   
  set(MPFR_FOUND TRUE)
  
else()  

  find_path(MPFR_INCLUDE_DIR 
            NAMES mpfr.h 
            PATHS ${CMAKE_SOURCE_DIR}/auxiliary/gmp/include
                  ENV MPFR_INC_DIR
  	        DOC "The directory containing the MPFR header files"
           )

  if ( MPFR_INCLUDE_DIR STREQUAL "${CMAKE_SOURCE_DIR}/auxiliary/gmp/include" )
    cache_set( MPFR_IN_CGAL_AUXILIARY TRUE )
  endif()
  
  if ( CGAL_AUTO_LINK_ENABLED )
  
    find_path(MPFR_LIBRARIES_DIR 
              NAMES "mpfr-${CGAL_TOOLSET}-mt.lib" "mpfr-${CGAL_TOOLSET}-mt-gd.lib"
              PATHS ${CMAKE_SOURCE_DIR}/auxiliary/gmp/lib
                    ENV MPFR_LIB_DIR
              DOC "Directory containing the MPFR library"
             ) 
    
  else()
  
    find_library(MPFR_LIBRARIES NAMES mpfr 
                 PATHS ENV MPFR_LIB_DIR
                 DOC "Path to the MPFR library"
                )
                
    if ( MPFR_LIBRARIES ) 
      get_filename_component(MPFR_LIBRARIES_DIR ${MPFR_LIBRARIES} PATH CACHE )
    endif()
    
  endif()  
  
  # Attempt to load a user-defined configuration for MPFR if couldn't be found
  if ( NOT MPFR_INCLUDE_DIR OR NOT MPFR_LIBRARIES_DIR )
    include( MPFRConfig OPTIONAL )
  endif()
  
  find_package_handle_standard_args(MPFR "DEFAULT_MSG" MPFR_INCLUDE_DIR MPFR_LIBRARIES_DIR)
  
endif()

