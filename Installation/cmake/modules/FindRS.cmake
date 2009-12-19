find_package( MPFI )

if( MPFI_FOUND )

  include( ${MPFI_USE_FILE} )

  find_path(RS_INCLUDE_DIR
            NAMES rs_exports.h
            PATHS ENV RS_INC_DIR
            DOC "The directory containing the RS include files"
           )

  find_path(RS3_INCLUDE_DIR
            NAMES rs3_fncts.h
            PATHS ENV RS_INC_DIR
            DOC "The directory containing the RS3 include files"
           )

  find_library(RS_LIBRARIES
               NAMES rsexport
               PATHS ENV RS_LIB_DIR
               DOC "Path to the RS library"
              )

  find_library(RS3_LIBRARIES
               NAMES rs3
               PATHS ENV RS_LIB_DIR
               DOC "Path to the RS3 library"
              )

  if( RS_INCLUDE_DIR AND RS_LIBRARIES )
    set(RS_FOUND TRUE)
  endif( RS_INCLUDE_DIR AND RS_LIBRARIES )

  if( RS3_INCLUDE_DIR AND RS3_LIBRARIES )
    set(RS3_FOUND TRUE)
  endif( RS3_INCLUDE_DIR AND RS3_LIBRARIES )

  if( RS_LIBRARIES )
    get_filename_component(RS_LIBRARIES_DIR ${RS_LIBRARIES} PATH CACHE )
  endif( RS_LIBRARIES )

  if( NOT RS_INCLUDE_DIR OR NOT RS_LIBRARIES_DIR )
    include( RSConfig OPTIONAL )
  endif( NOT RS_INCLUDE_DIR OR NOT RS_LIBRARIES_DIR )

  include(CGAL_FindPackageHandleStandardArgs)

  find_package_handle_standard_args( RS
                                     "DEFAULT_MSG"
                                     RS_LIBRARIES
                                     RS_INCLUDE_DIR )
else( MPFI_FOUND )

  message( STATUS "RS requires MPFI" )
  set( RS_FOUND FALSE )

endif( MPFI_FOUND )

if(RS_FOUND)
  set(RS_USE_FILE "CGAL_UseRS")
endif(RS_FOUND)
