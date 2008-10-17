if (TAUCS_INCLUDE_DIR AND TAUCS_LIBRARIES_DIR ) 
   
  set(TAUCS_FOUND TRUE)
  
else()  
 
  include(CGAL_Locate_CGAL_TAUCS)
  if( CGAL_TAUCS_FOUND )
  
     set( TAUCS_INCLUDE_DIR   "${CGAL_TAUCS_INCLUDE_DIR}"   CACHE FILEPATH "Include directories for the TAUCS libraries" )
     
     set( TAUCS_LIBRARIES_DIR "${CGAL_TAUCS_LIBRARIES_DIR}" CACHE FILEPATH "Lib directories for the TAUCS libraries")
     
     set( TAUCS_LIBRARIES "${CGAL_TAUCS_LIBRARIES_DIR}/libtaucs.a;${CGAL_TAUCS_DIR}/external/lib/${CGAL_TAUCS_PLATFORM}/libmetis.a;${CGAL_TAUCS_DIR}/external/lib/${CGAL_TAUCS_PLATFORM}/libatlas.a;${CGAL_TAUCS_DIR}/external/lib/${CGAL_TAUCS_PLATFORM}/libg2c.so"
          CACHE FILEPATH "The TAUCS libraries"                    
        )
     
  else()
  
    find_path(TAUCS_INCLUDE_DIR 
              NAMES taucs.h 
              PATHS ENV TAUCS_INC_DIR
              DOC "The directory containing the TAUCS header files"
           )
           
    find_library(TAUCS_LIBRARIES NAMES "taucs"
                 PATHS ENV TAUCS_LIB_DIR
                 DOC "Path to the TAUCS library"
                )
                
    if ( TAUCS_LIBRARIES ) 
      get_filename_component(TAUCS_LIBRARIES_DIR ${TAUCS_LIBRARIES} PATH CACHE)
    endif()
    
  endif()
  
  if ( TAUCS_INCLUDE_DIR AND TAUCS_LIBRARIES_DIR)
    set(TAUCS_FOUND TRUE)
  endif()
  
endif()
