find_package( GMP QUIET )
find_package( MPFR QUIET )

if( GMP_FOUND AND MPFR_FOUND )

  if( MPFI_INCLUDE_DIR AND MPFI_LIBRARIES )
    set( MPFI_FOUND TRUE )
  endif( MPFI_INCLUDE_DIR AND MPFI_LIBRARIES )

  find_path(MPFI_INCLUDE_DIR NAMES mpfi.h
            HINTS ENV MPFI_INC_DIR
                  ENV MPFI_DIR
            PATHS ${GMP_INCLUDE_DIR_SEARCH}
            PATH_SUFFIXES include
            DOC "The directory containing the MPFI header files"
           )

  find_library(MPFI_LIBRARIES NAMES mpfi
               HINTS ENV MPFI_LIB_DIR
                     ENV MPFI_DIR
               PATHS ${GMP_LIBRARIES_DIR_SEARCH}
               PATH_SUFFIXES lib
               DOC "Directory containing the MPFI library"
               )

  if( MPFI_LIBRARIES )
    get_filename_component(MPFI_LIBRARIES_DIR ${MPFI_LIBRARIES} PATH CACHE )
  endif( MPFI_LIBRARIES )

  if( NOT MPFI_INCLUDE_DIR OR NOT MPFI_LIBRARIES_DIR )
    include( MPFIConfig OPTIONAL )
  endif( NOT MPFI_INCLUDE_DIR OR NOT MPFI_LIBRARIES_DIR )

  include(FindPackageHandleStandardArgs)

  find_package_handle_standard_args( MPFI
                                     "DEFAULT_MSG"
                                     MPFI_LIBRARIES
                                     MPFI_INCLUDE_DIR )

else( GMP_FOUND AND MPFR_FOUND )

  message( STATUS "MPFI needs GMP and MPFR" )

endif( GMP_FOUND AND MPFR_FOUND )

if( MPFI_FOUND )
  get_dependency_version( MPFR )
  IS_VERSION_LESS("${MPFR_VERSION}" "4.0.0" _MPFR_OLD)

  get_dependency_version( MPFI )
  IS_VERSION_LESS("${MPFI_VERSION}" "1.5.2" _MPFI_OLD)

  if( ( _MPFR_OLD AND NOT _MPFI_OLD ) OR ( NOT _MPFR_OLD AND _MPFI_OLD ) )

    message(
      STATUS
      "MPFI<1.5.2 requires MPFR<4.0.0; MPFI>=1.5.2 requires MPFR>=4.0.0" )

    set( MPFI_FOUND FALSE )

  else( ( _MPFR_OLD AND NOT _MPFI_OLD ) OR ( NOT _MPFR_OLD AND _MPFI_OLD ) )

    set( MPFI_USE_FILE "CGAL_UseMPFI" )

  endif( ( _MPFR_OLD AND NOT _MPFI_OLD ) OR ( NOT _MPFR_OLD AND _MPFI_OLD ) )

endif( MPFI_FOUND )
