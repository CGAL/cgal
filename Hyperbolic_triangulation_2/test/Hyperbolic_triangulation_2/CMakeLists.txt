# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.15)
project( Hyperbolic_triangulation_2_Tests )


find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

          include( ${CGAL_USE_FILE} )

          include( CGAL_CreateSingleSourceCGALProgram )

          create_single_source_cgal_program( "ht2_test_clear.cpp" )
          target_link_libraries(ht2_test_clear PRIVATE CGAL::CGAL_Core)
          create_single_source_cgal_program( "ht2_test_locate.cpp" )
          target_link_libraries(ht2_test_locate PRIVATE CGAL::CGAL_Core)
          create_single_source_cgal_program( "ht2_test_remove.cpp" )
          target_link_libraries(ht2_test_remove PRIVATE CGAL::CGAL_Core)
          create_single_source_cgal_program( "ht2_test_swap.cpp" )
          target_link_libraries(ht2_test_swap PRIVATE CGAL::CGAL_Core)
          create_single_source_cgal_program( "ht2_test_copy.cpp" )
          target_link_libraries(ht2_test_copy PRIVATE CGAL::CGAL_Core)
          create_single_source_cgal_program( "ht2_test_hyperbolic_circulator.cpp" )
          target_link_libraries(ht2_test_hyperbolic_circulator PRIVATE CGAL::CGAL_Core)
          create_single_source_cgal_program( "ht2_test_insert_degenerate.cpp" )
          target_link_libraries(ht2_test_insert_degenerate PRIVATE CGAL::CGAL_Core)

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

