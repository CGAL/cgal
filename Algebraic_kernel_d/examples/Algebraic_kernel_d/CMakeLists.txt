cmake_minimum_required(VERSION 3.1...3.15)
project( Algebraic_kernel_d_Examples )



find_package(CGAL QUIET COMPONENTS Core)

if ( CGAL_FOUND )
  find_package(MPFI QUIET)
endif()

if( CGAL_FOUND AND MPFI_FOUND)

  include( ${CGAL_USE_FILE} )
  include( ${MPFI_USE_FILE} )
  include( CGAL_VersionUtils )

  create_single_source_cgal_program( "Compare_1.cpp" )
  target_link_libraries(Compare_1 PRIVATE CGAL::CGAL_Core)
  create_single_source_cgal_program( "Construct_algebraic_real_1.cpp" )
  target_link_libraries(Construct_algebraic_real_1 PRIVATE CGAL::CGAL_Core)
  create_single_source_cgal_program( "Isolate_1.cpp" )
  target_link_libraries(Isolate_1 PRIVATE CGAL::CGAL_Core)
  create_single_source_cgal_program( "Sign_at_1.cpp" )
  target_link_libraries(Sign_at_1 PRIVATE CGAL::CGAL_Core)
  create_single_source_cgal_program( "Solve_1.cpp" )
  target_link_libraries(Solve_1 PRIVATE CGAL::CGAL_Core)

else()

    message(STATUS "This program requires the CGAL library and MPFI, and will not be compiled.")

endif()
