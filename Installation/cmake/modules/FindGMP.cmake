# Try to find the GMP libraries
# GMP_FOUND - system has GMP lib
# GMP_INCLUDE_DIR - the GMP include directory
# GMP_LIBRARIES_DIR - Directory where the GMP libraries are located
# GMP_LIBRARIES - the GMP libraries
# GMP_IN_CGAL_AUXILIARY - TRUE if the GMP found is the one distributed with CGAL in the auxiliary folder

# TODO: support MacOSX

include(CGAL_FindPackageHandleStandardArgs)
include(CGAL_GeneratorSpecificSettings)

# Is it already configured?
if (GMP_INCLUDE_DIR AND GMP_LIBRARIES_DIR ) 
   
  set(GMP_FOUND TRUE)
  
else()  

  find_path(GMP_INCLUDE_DIR 
            NAMES gmp.h 
            PATHS ${CMAKE_SOURCE_DIR}/auxiliary/gmp/include
                  ENV GMP_INC_DIR
  	        DOC "The directory containing the GMP header files"
           )

  if ( GMP_INCLUDE_DIR STREQUAL "${CMAKE_SOURCE_DIR}/auxiliary/gmp/include" )
    cache_set( GMP_IN_CGAL_AUXILIARY TRUE )
  endif()
  
  if ( CGAL_AUTO_LINK_ENABLED )
  
    find_path(GMP_LIBRARIES_DIR 
              NAMES "gmp-${CGAL_TOOLSET}-mt.lib" "gmp-${CGAL_TOOLSET}-mt-gd.lib"
              PATHS ${CMAKE_SOURCE_DIR}/auxiliary/gmp/lib
                    ENV GMP_LIB_DIR
              DOC "Directory containing the GMP library"
             ) 
    
  else()
  
    find_library(GMP_LIBRARIES NAMES gmp 
                 PATHS ENV GMP_LIB_DIR
                 DOC "Path to the GMP library"
                )
                
    if ( GMP_LIBRARIES ) 
      get_filename_component(GMP_LIBRARIES_DIR ${GMP_LIBRARIES} PATH CACHE )
    endif()
    
  endif()  
    
  # Attempt to load a user-defined configuration for GMP if couldn't be found
  if ( NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARIES_DIR )
    include( GMPConfig OPTIONAL )
  endif()
  
  find_package_handle_standard_args(GMP "DEFAULT_MSG" GMP_INCLUDE_DIR GMP_LIBRARIES_DIR)
  
endif()
