cmake_minimum_required(VERSION 3.1...3.15)
project( FastEnvelope_Examples )


include_directories( BEFORE  "c:/3rdPartyLibs/fast-envelope-VC/include" )

link_directories ( "c:/3rdPartyLibs/fast-envelope-VC/lib"  "c:/3rdPartyLibs/fast-envelope-VC"  "c:/3rdPartyLibs/fast-envelope-VC/tbb"   )

find_package(CGAL QUIET)

if ( CGAL_FOUND )

find_package(Eigen3 3.2.0) #(requires 3.2.0 or greater)
include(CGAL_Eigen_support)

#  add_definitions(-DQT_NO_KEYWORDS)

    create_single_source_cgal_program( "aabb.cpp" )
    create_single_source_cgal_program( "fast.cpp" )
    create_single_source_cgal_program( "fastE.cpp" )

    target_link_libraries( fast PUBLIC CGAL::Eigen_support)

    target_link_libraries( fastE PUBLIC CGAL::Eigen_support)
  target_link_libraries( fastE PUBLIC ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} FastEnvelope IndirectPredicates geogram)
else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
