cmake_minimum_required(VERSION 3.12...3.31)
project(Arrangement_on_curve_1_Examples)
find_package(CGAL OPTIONAL_COMPONENTS Qt6)

# create a target per cppfile
file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})
  create_single_source_cgal_program("${cppfile}")
endforeach()

if(CGAL_Qt6_FOUND)
  target_link_libraries(arr_on_circle_segment PRIVATE CGAL::CGAL_Basic_viewer)
else()
  message(
    STATUS
      "NOTICE: Several examples require Qt and will not be compiled."
  )
endif()
