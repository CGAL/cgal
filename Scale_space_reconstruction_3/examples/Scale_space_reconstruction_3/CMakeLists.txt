project( Scale_space_reconstruction_3_Examples )

cmake_minimum_required(VERSION 2.8.11)

find_package( CGAL QUIET)

if ( CGAL_FOUND )
  include( ${CGAL_USE_FILE} )

  # Activate concurrency?
  option(ACTIVATE_CONCURRENCY
         "Enable concurrency"
         ON)
 
  if( ACTIVATE_CONCURRENCY )
    find_package( TBB )
    if( NOT TBB_FOUND )
      message( STATUS "NOTICE: Intel TBB NOT found! The example is faster if TBB is linked." )
    endif()
  endif()

  include( CGAL_CreateSingleSourceCGALProgram )

  find_package( Eigen3 3.1.0 )
  if( EIGEN3_FOUND )
    include( ${EIGEN3_USE_FILE} )
    create_single_source_cgal_program( "scale_space.cpp" )
    create_single_source_cgal_program( "scale_space_incremental.cpp" )
    create_single_source_cgal_program( "scale_space_manifold.cpp" )
    create_single_source_cgal_program( "scale_space_advancing_front.cpp" )
    if(ACTIVATE_CONCURRENCY AND TBB_FOUND)
      CGAL_target_use_TBB(scale_space)
      CGAL_target_use_TBB(scale_space_incremental)
      CGAL_target_use_TBB(scale_space_manifold)
      CGAL_target_use_TBB(scale_space_advancing_front)
    endif()
  else()
    message( STATUS "NOTICE: The example needs Eigen 3.1 (or greater) and will not be compiled." )
  endif()
else()
  
  message( STATUS "This program requires the CGAL library, and will not be compiled." )
  
endif()
