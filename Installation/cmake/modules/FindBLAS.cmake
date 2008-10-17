if ( BLAS_LIBRARIES ) 
   
  set(BLAS_FOUND TRUE)
  
else()  
 
  include(CGAL_Locate_CGAL_TAUCS)
  if( CGAL_TAUCS_FOUND )
                                                                                                                                                    
     set( BLAS_LIBRARIES    "${CGAL_TAUCS_DIR}/external/lib/${CGAL_TAUCS_PLATFORM}/libcblas.a;${CGAL_TAUCS_DIR}/external/lib/${CGAL_TAUCS_PLATFORM}/libf77blas.a" CACHE FILEPATH "The BLAS libraries" )
     set( BLAS_DEFINITIONS "-DCGAL_USE_F2C -DCGAL_USE_CBLASWRAP"                                                                                                  CACHE STRING   "Definitions for the BLAS" )
     
  else()
  
    set( BLAS_LIBRARIES   "$ENV{BLAS_LIBRARIES}"   CACHE FILEPATH "The BLAS libraries" )
    set( BLAS_DEFINITIONS "$ENV{BLAS_DEFINITIONS}" CACHE STRING   "Definitions for the BLAS" )

  endif()
  
  if ( BLAS_LIBRARIES )
    set(BLAS_FOUND TRUE)
  endif()
  
endif()
