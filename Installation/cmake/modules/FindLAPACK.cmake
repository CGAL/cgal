if ( LAPACK_LIBRARIES ) 
   
  set(LAPACK_FOUND TRUE)
  
else()  
 
  include(CGAL_Locate_CGAL_TAUCS)
  if( CGAL_TAUCS_FOUND )
  
     set( LAPACK_LIBRARIES "${CGAL_TAUCS_DIR}/external/lib/${CGAL_TAUCS_PLATFORM}/liblapack.a" CACHE FILEPATH "The LAPACK libraries" )
     
  else()
  
    set( LAPACK_LIBRARIES   $ENV{LAPACK_LIBRARIES}   CACHE FILEPATH "The LAPACK libraries" )
    set( LAPACK_DEFINITIONS $ENV{LAPACK_DEFINITIONS} CACHE STRING   "Definitions for the LAPACK libraries" )
              
  endif()
  
  if ( LAPACK_LIBRARIES )
    set(LAPACK_FOUND TRUE)
  endif()
  
endif()
