# Try to find the GMP libraries
# GMP_FOUND - system has GMP lib
# GMP_INCLUDE_DIR - the GMP include directory
# GMP_LIBRARIES_DIR - Directory where the GMP libraries are located
# GMP_LIBRARIES - the GMP libraries
# GMP_IN_CGAL_AUXILIARY - TRUE if the GMP found is the one distributed with CGAL in the auxiliary folder

# TODO: support MacOSX

include(FindPackageHandleStandardArgs)
include(GeneratorSpecificSettings)
include(ReadLines)
include(FindMatchingItem)

# Is it already configured?
if (GMP_INCLUDE_DIR AND GMP_LIBRARIES_DIR ) 
   
  set(GMP_FOUND TRUE)
  
else()  

  # Look first for the GMP distributed with CGAL in auxiliary/gmp
  find_path(GMP_INCLUDE_DIR 
            NAMES gmp.h 
            PATHS ${CGAL_SOURCE_DIR}/auxiliary/gmp/include
                  ENV GMP_INC_DIR
  	        DOC "The directory containing the GMP header files"
           )

  if ( GMP_INCLUDE_DIR STREQUAL "${CGAL_SOURCE_DIR}/auxiliary/gmp/include" )
    set( GMP_IN_CGAL_AUXILIARY TRUE CACHE INTERNAL "" )
  endif()
  
  if ( AUTO_LINK_ENABLED )
  
    find_path(GMP_LIBRARIES_DIR 
              NAMES "gmp-${TOOLSET}-mt.lib" "gmp-${TOOLSET}-mt-gd.lib"
              PATHS ${CGAL_SOURCE_DIR}/auxiliary/gmp/lib
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
  
  if ( GMP_FOUND )
  
    readlines(${GMP_INCLUDE_DIR}/gmp.h GMP_H_FILE)
  
    find_matching_item(GMP_H_FILE "__GNU_MP_VERSION "            __GNU_MP_VERSION_LINE )
    find_matching_item(GMP_H_FILE "__GNU_MP_VERSION_MINOR "      __GNU_MP_VERSION_MINOR_LINE )
    find_matching_item(GMP_H_FILE "__GNU_MP_VERSION_PATCHLEVEL " __GNU_MP_VERSION_PATCHLEVEL_LINE )
    
    string( REPLACE " " ";" __GNU_MP_VERSION_LINE_LIST            ${__GNU_MP_VERSION_LINE}            )
    string( REPLACE " " ";" __GNU_MP_VERSION_MINOR_LINE_LIST      ${__GNU_MP_VERSION_MINOR_LINE}      )
    string( REPLACE " " ";" __GNU_MP_VERSION_PATCHLEVEL_LINE_LIST ${__GNU_MP_VERSION_PATCHLEVEL_LINE} )
    
    list( GET __GNU_MP_VERSION_LINE_LIST            2 __GNU_MP_VERSION )
    list( GET __GNU_MP_VERSION_MINOR_LINE_LIST      2 __GNU_MP_VERSION_MINOR )
    list( GET __GNU_MP_VERSION_PATCHLEVEL_LINE_LIST 2 __GNU_MP_VERSION_PATCHLEVEL )
    
    set( GMP_VERSION "${__GNU_MP_VERSION}.${__GNU_MP_VERSION_MINOR}.${__GNU_MP_VERSION_PATCHLEVEL}" )

    message( STATUS "USING GMP_VERSION ${GMP_VERSION}" )
    
  endif()
  
endif()

#mark_as_advanced(GMP_INCLUDE_DIR)
#mark_as_advanced(GMP_LIBRARIES)
#mark_as_advanced(GMP_LIBRARIES_DIR)
