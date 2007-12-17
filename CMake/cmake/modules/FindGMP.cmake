# Try to find the GMP libraries
# GMP_FOUND - system has GMP lib
# GMP_INCLUDE_DIR - the GMP include directory
# GMP_LIBRARIES - Libraries needed to use GMP

# TODO: support MacOSX

include(FindPackageHandleStandardArgs)
include(GeneratorSpecificSettings)

if (GMP_INCLUDE_DIR AND GMP_LIBRARIES)
  # Already in cache, be silent
  set(GMP_FIND_QUIETLY TRUE)
endif()

set(GMP_INCLUDE_DIR_SEARCH   ${CGAL_SOURCE_DIR}/auxiliary/gmp/include)
set(GMP_LIBRARIES_DIR_SEARCH ${CGAL_SOURCE_DIR}/auxiliary/gmp/lib)

find_path(GMP_INCLUDE_DIR NAMES gmp.h PATHS
	        ${GMP_INCLUDE_DIR_SEARCH}
	        DOC "The directory containing the GMP include files"
         )

if ( AUTO_LINK_ENABLED )
    
  set(GMP_LIBRARIES "" )
  find_path(GMP_LIBRARIES_DIR 
            NAMES "gmp${TOOLSET}-mt.lib" "gmp${TOOLSET}-mt-gd.lib" "gmp${TOOLSET}-mt-o.lib" "gmp${TOOLSET}-mt-g.lib"
            PATHS ${GMP_LIBRARIES_DIR_SEARCH}
            DOC "Directory containing the GMP library"
           ) 
    
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(GMP "DEFAULT_MSG" GMP_INCLUDE_DIR GMP_LIBRARIES_DIR)
  
else()

  find_library(GMP_LIBRARIES NAMES gmp 
               PATHS ${GMP_LIBRARIES_DIR_SEARCH}
  	           DOC "Path to the GMP library"
              )
              
  get_filename_component(GMP_LIBRARIES_DIR ${GMP_LIBRARIES} PATH)
  
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(GMP "DEFAULT_MSG" GMP_INCLUDE_DIR GMP_LIBRARIES )
  
endif()

