# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Voronoi_diagram_2_Examples )

cmake_minimum_required(VERSION 3.1)

find_package(CGAL COMPONENTS Qt5)

if(CGAL_Qt5_FOUND)
  add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
endif()

if ( CGAL_FOUND )

  # create a target per cppfile
  file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
  foreach(cppfile ${cppfiles})
    create_single_source_cgal_program( "${cppfile}" )
  endforeach()

  if(CGAL_Qt5_FOUND )
    target_link_libraries(vd_2_point_location PUBLIC CGAL::CGAL_Qt5)
    target_link_libraries(vd_2_point_location_sdg_linf PUBLIC CGAL::CGAL_Qt5)
    target_link_libraries(draw_voronoi_diagram_2 PUBLIC CGAL::CGAL_Qt5)
  endif()


else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

