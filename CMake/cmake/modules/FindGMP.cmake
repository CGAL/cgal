# Try to find the GMP libraries
# GMP_FOUND - system has GMP lib
# GMP_INCLUDE_DIR - the GMP include directory
# GMP_LIBRARIES_DIR - Directory where the GMP libraries are located
# GMP_LIBRARIES - the GMP libraries
# GMP_IN_AUXILIARY - TRUE if the GMP found is the one distributed with CGAL in the auxiliary folder

# TODO: support MacOSX

include(FindPackageHandleStandardArgs)
include(GeneratorSpecificSettings)

# Is it already configured?
if (GMP_INCLUDE_DIR AND GMP_LIBRARIES ) 
   
  set(GMP_FOUND TRUE)
  
else()  

  # Look first for the GMP distributed with CGAL in auxiliary/gmp
  find_path(GMP_INCLUDE_DIR 
            NAMES gmp.h 
            PATHS ${CGAL_SOURCE_DIR}/auxiliary/gmp/include
                  ENV GMP_INC_DIR
  	        DOC "The directory containing the GMP header files"
           )

  if ( GMP_INCLUDE_DIR ) 
     if ( GMP_INCLUDE_DIR STREQUAL "${CGAL_SOURCE_DIR}/auxiliary/gmp/include" )
       set( GMP_IN_CGAL_AUXILIARY TRUE )
       set( GMP_LIB_SEARCH_PATHS ${CGAL_SOURCE_DIR}/auxiliary/gmp/lib )
     endif()
     
    if ( AUTO_LINK_ENABLED )
      set( GMP_NAMES "gmp${TOOLSET}-mt" "gmp${TOOLSET}-mt-gd" "gmp${TOOLSET}-mt-o" "gmp${TOOLSET}-mt-g" )
    else()
      set( GMP_NAMES "gmp" )
    endif()  
    
    find_library(GMP_LIBRARIES 
                 NAMES ${GMP_NAMES} 
                 PATHS ${GMP_LIB_SEARCH_PATHS}
                       ENV GMP_LIB_DIR
                 DOC "Path to the GMP library"
                )
     
    if ( GMP_LIBRARIES ) 
      get_filename_component(GMP_LIBRARIES_DIR ${GMP_LIBRARIES} PATH)
    endif()
    
  endif()  
  
  # Attempt to load a user-defined configuration for GMP if couldn't be found
  if ( NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARIES )
    include( GMPConfig OPTIONAL )
  endif()
  
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(GMP "DEFAULT_MSG" GMP_INCLUDE_DIR GMP_LIBRARIES)
  
endif()


