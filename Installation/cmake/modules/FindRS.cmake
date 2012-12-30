# RS needs GMP 4.2 or newer, this script will fail if an old version is
# detected

if( NOT GMP_FOUND )
  find_package( GMP )
endif()

if( NOT MPFI_FOUND )
  find_package( MPFI )
endif()

if( MPFI_FOUND )

  include( CGAL_VersionUtils )

  find_path(RS_INCLUDE_DIR
            NAMES rs_exports.h
            HINTS ENV RS_INC_DIR
                  ENV RS_DIR
            PATH_SUFFIXES include
            DOC "The directory containing the RS include files"
           )

  find_library(RS_LIBRARIES
               NAMES rsexport_rs
               HINTS ENV RS_LIB_DIR
                     ENV RS_DIR
               PATH_SUFFIXES lib
               DOC "Path to the RS library"
              )

  get_dependency_version( GMP )

  IS_VERSION_LESS("${GMP_VERSION}" "4.2.0" _IS_GMP_VERSION_TO_LOW)

  if(_IS_GMP_VERSION_TO_LOW)

    message( STATUS
      "RS needs GMP>=4.2. Your GMP version is ${GMP_VERSION}." )

  else(_IS_GMP_VERSION_TO_LOW)

    if( RS_INCLUDE_DIR AND RS_LIBRARIES )
      set(RS_FOUND TRUE)
    endif( RS_INCLUDE_DIR AND RS_LIBRARIES )

    if( RS_LIBRARIES )
      get_filename_component(RS_LIBRARIES_DIR ${RS_LIBRARIES} PATH CACHE )
    endif( RS_LIBRARIES )

    include(FindPackageHandleStandardArgs)

    find_package_handle_standard_args( RS
                                       "DEFAULT_MSG"
                                       RS_LIBRARIES
                                       RS_INCLUDE_DIR )

  endif(_IS_GMP_VERSION_TO_LOW)

else( MPFI_FOUND )

  message( STATUS "RS requires MPFI" )
  set( RS_FOUND FALSE )

endif( MPFI_FOUND )

if(RS_FOUND)
  set(RS_USE_FILE "CGAL_UseRS")
endif(RS_FOUND)
