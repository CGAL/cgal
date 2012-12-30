# RS needs GMP 4.2 or newer, this script will fail if an old version is
# detected

find_package( GMP )
find_package( MPFI )
find_package( RS )

if ( RS_FOUND )

if( MPFI_FOUND )

  include( CGAL_VersionUtils )

  find_path(RS3_INCLUDE_DIR
            NAMES rs3_fncts.h
            HINTS ENV RS3_INC_DIR
                  ENV RS3_DIR
                  ENV RS_INC_DIR
                  ENV RS_DIR
            PATHS ${RS_INCLUDE_DIR}
            PATH_SUFFIXES include
            DOC "The directory containing the RS3 include files"
           )

  find_library(RS3_LIBRARIES
               NAMES rs3
               HINTS ENV RS3_LIB_DIR
                     ENV RS3_DIR
                     ENV RS_LIB_DIR
                     ENV RS_DIR
               PATHS ${RS_LIBRARIES_DIR}
               PATH_SUFFIXES lib
               DOC "Path to the RS3 library"
              )

  get_dependency_version( GMP )

  IS_VERSION_LESS("${GMP_VERSION}" "4.2.0" _IS_GMP_VERSION_TO_LOW)

  if(_IS_GMP_VERSION_TO_LOW)

    message( STATUS
      "RS3 needs GMP>=4.2. Your GMP version is ${GMP_VERSION}." )

  else(_IS_GMP_VERSION_TO_LOW)

    if( RS3_INCLUDE_DIR AND RS3_LIBRARIES )
      set(RS3_FOUND TRUE)
    endif( RS3_INCLUDE_DIR AND RS3_LIBRARIES )

    if( RS3_LIBRARIES )
      get_filename_component(RS3_LIBRARIES_DIR ${RS3_LIBRARIES} PATH CACHE )
    endif( RS3_LIBRARIES )

    if( NOT RS3_INCLUDE_DIR OR NOT RS3_LIBRARIES_DIR )
      include( RS3Config OPTIONAL )
    endif( NOT RS3_INCLUDE_DIR OR NOT RS3_LIBRARIES_DIR )

    include(FindPackageHandleStandardArgs)

    find_package_handle_standard_args( RS3
                                       "DEFAULT_MSG"
                                       RS3_LIBRARIES
                                       RS3_INCLUDE_DIR )

  endif(_IS_GMP_VERSION_TO_LOW)

else( MPFI_FOUND )

  message( STATUS "RS3 requires MPFI" )
  set( RS3_FOUND FALSE )

endif( MPFI_FOUND )

else( RS_FOUND )

  message( STATUS "RS3 requires RS" )
  set( RS3_FOUND FALSE )

endif( RS_FOUND )

if(RS3_FOUND)
  set(RS3_USE_FILE "CGAL_UseRS3")
endif(RS3_FOUND)
