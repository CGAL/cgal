find_package( GMP QUIET )

if( GMP_FOUND )
  
  if( MPFI_INCLUDE_DIR AND MPFI_LIBRARIES )
    set( MPFI_FOUND TRUE )
  endif( MPFI_INCLUDE_DIR AND MPFI_LIBRARIES )

  find_path(MPFI_INCLUDE_DIR NAMES mpfi.h
            HINTS
            $ENV{MPFI_INC_DIR}
            PATHS 
            ${GMP_INCLUDE_DIR_SEARCH}
            DOC "The directory containing the MPFI header files"
           )

  find_library(MPFI_LIBRARIES NAMES mpfi
               HINTS
               $ENV{MPFI_LIB_DIR}
               PATHS 
               ${GMP_LIBRARIES_DIR_SEARCH}
               DOC "Directory containing the MPFI library"
               )

  if( MPFI_LIBRARIES )
    get_filename_component(MPFI_LIBRARIES_DIR ${MPFI_LIBRARIES} PATH CACHE )
  endif( MPFI_LIBRARIES )

  if( NOT MPFI_INCLUDE_DIR OR NOT MPFI_LIBRARIES_DIR )
    include( MPFIConfig OPTIONAL )
  endif( NOT MPFI_INCLUDE_DIR OR NOT MPFI_LIBRARIES_DIR )

  include(FindPackageHandleStandardArgs)

  find_package_handle_standard_args(MPFI "DEFAULT_MSG" MPFI_LIBRARIES MPFI_INCLUDE_DIR )

else( GMP_FOUND )

  message( STATUS "MPFI needs GMP and MPFR" )

endif( GMP_FOUND )

if( MPFI_FOUND )
  set( MPFI_USE_FILE "CGAL_UseMPFI" )
endif( MPFI_FOUND )
