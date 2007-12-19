# Try to find the MPFR libraries
# MPFR_FOUND - system has MPFR lib
# MPFR_INCLUDE_DIR - the MPFR include directory
# MPFR_LIBRARIES - Libraries needed to use MPFR

# TODO: support Windows and MacOSX

include(FindPackageHandleStandardArgs)

# MPFR needs GMP
find_package(GMP QUIET)
if(GMP_FOUND)

  if (MPFR_INCLUDE_DIR AND MPFR_LIBRARIES)
    # Already in cache, be silent
    set(MPFR_FIND_QUIETLY TRUE)
  endif()

  find_path(MPFR_INCLUDE_DIR NAMES mpfr.h
           PATHS ${GMP_INCLUDE_DIR_SEARCH}
           DOC "The directory containing the MPFR include files"
           )

  if ( AUTO_LINK_ENABLED )
  
    set(MPFR_LIBRARIES "" )
    
    find_path(MPFR_LIBRARIES_DIR 
              NAMES "mpfr${TOOLSET}-mt.lib" "mpfr${TOOLSET}-mt-gd.lib" "mpfr${TOOLSET}-mt-o.lib" "mpfr${TOOLSET}-mt-g.lib"
              PATHS ${GMP_LIBRARIES_DIR_SEARCH}
              DOC "Directory containing the MPFR library"
             ) 
             
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(MPFR "DEFAULT_MSG" MPFR_LIBRARIES_DIR MPFR_INCLUDE_DIR )
    
  else()
  
    find_library(MPFR_LIBRARIES NAMES mpfr 
                 PATHS ${GMP_LIBRARIES_DIR_SEARCH}
                 DOC "Path to the MPFR library"
                )
                
    get_filename_component(MPFR_LIBRARIES_DIR ${MPFR_LIBRARIES} PATH)
    
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(MPFR "DEFAULT_MSG" MPFR_INCLUDE_DIR MPFR_LIBRARIES )
    
  endif()

endif(GMP_FOUND)
