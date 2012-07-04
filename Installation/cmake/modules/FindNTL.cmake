# NTL needs GMP 3.1.1 or newer, this script will fail if an old version is
# detected

find_package( GMP REQUIRED )

if( (TARGET CGAL AND NOT WITH_GMP) OR NOT GMP_FOUND )

  message( FATAL_ERROR "NTL requires GMP" )
  set( NTL_FOUND FALSE )

else( (TARGET CGAL AND NOT WITH_GMP) OR NOT GMP_FOUND )

  include( CGAL_VersionUtils )

  get_dependency_version( GMP )

  IS_VERSION_LESS("${GMP_VERSION}" "3.1.1" _IS_GMP_VERSION_TOO_LOW)

  if( _IS_GMP_VERSION_TOO_LOW )

    message( FATAL_ERROR, "NTL needs GMP>=3.1.1. Your GMP version is ${GMP_VERSION}." )

  else( _IS_GMP_VERSION_TOO_LOW )

    find_path(NTL_INCLUDE_DIR
              NAMES NTL/ZZ.h
              HINTS
              $ENV{NTL_INC_DIR}
              DOC "The directory containing the NTL include files"
             )

    find_library(NTL_LIBRARY
                 NAMES ntl
                 HINTS
                 $ENV{NTL_LIB_DIR}
                 DOC "Path to the NTL library"
                )

    if ( NTL_INCLUDE_DIR AND NTL_LIBRARY ) 
      
       #check version

       set( NTL_VERSION_H "${NTL_INCLUDE_DIR}/NTL/version.h" )

       if ( EXISTS ${NTL_VERSION_H} )

         file( READ "${NTL_VERSION_H}" NTL_VERSION_H_CONTENTS )
  
         string( REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" CGAL_NTL_VERSION "${NTL_VERSION_H_CONTENTS}" )
  
         #message( STATUS "DETECTED NTL_VERSION = '${CGAL_NTL_VERSION}'" )

         IS_VERSION_GREATER( "${CGAL_NTL_VERSION}" "5.0.0" _IS_NTL_VERSION_GREATER )

         if ( _IS_NTL_VERSION_GREATER )
           set( NTL_FOUND TRUE )
         else ( _IS_NTL_VERSION_GREATER )
           set( NTL_FOUND FALSE )
         endif ( _IS_NTL_VERSION_GREATER )


       endif (EXISTS ${NTL_VERSION_H} )

    endif ( NTL_INCLUDE_DIR AND NTL_LIBRARY ) 

    if ( NTL_FOUND )

      #message( STATUS "Found NTL in version '${CGAL_NTL_VERSION}'" )
      set ( NTL_INCLUDE_DIRS ${NTL_INCLUDE_DIR} )
      set ( NTL_LIBRARIES ${NTL_LIBRARY} )

      get_filename_component(NTL_LIBRARIES_DIR ${NTL_LIBRARIES} PATH CACHE )

      include(FindPackageHandleStandardArgs)

      find_package_handle_standard_args( NTL
                                         DEFAULT_MSG
                                         NTL_LIBRARY
                                         NTL_INCLUDE_DIR )
    
      mark_as_advanced( NTL_INCLUDE_DIR NTL_LIBRARY )

      # TODO add flag to CGAL Polynomials

    endif( NTL_FOUND )


  endif( _IS_GMP_VERSION_TOO_LOW )

endif( (TARGET CGAL AND NOT WITH_GMP) OR NOT GMP_FOUND )

if ( NTL_FOUND )
#  if ( NOT NTL_FIND_QUIETLY )
#    message(STATUS "Found NTL: ${NTL_LIBRARY}")
#  endif (NOT NTL_FIND_QUIETLY )
else ( NTL_FOUND )
  if ( NTL_FIND_REQUIRED )
    message( FATAL_ERROR "Could not find NTL" )
  endif ( NTL_FIND_REQUIRED )
endif ( NTL_FOUND )
