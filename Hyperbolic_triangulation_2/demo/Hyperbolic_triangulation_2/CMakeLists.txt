cmake_minimum_required(VERSION 3.12...3.31)
project(Hyperbolic_triangulation_2_Demo)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Core Qt6)
find_package(LEDA QUIET)

# Find Qt6 itself
find_package(Qt6 QUIET COMPONENTS Widgets)

if(CGAL_Qt6_FOUND
   AND Qt6_FOUND
   AND (CGAL_Core_FOUND OR LEDA_FOUND))

   set(CMAKE_AUTOMOC ON)
   set(CMAKE_AUTOUIC ON)
   set(CMAKE_AUTORCC ON)

  qt_add_executable ( HDT2 HDT2.cpp HDT2.ui resources/Delaunay_triangulation_2.qrc )
  target_include_directories(HDT2 PRIVATE ./ ./include)
  add_to_cached_list( CGAL_EXECUTABLE_TARGETS HDT2 )
  target_link_libraries ( HDT2 PRIVATE CGAL::CGAL CGAL::CGAL_Qt6 Qt6::Widgets)
  if(CGAL_Core_FOUND)
    target_link_libraries ( HDT2 PRIVATE CGAL::CGAL_Core)
  else()
    target_link_libraries ( HDT2 PRIVATE ${LEDA_LIBRARIES})
  endif()

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test( HDT2 )
else()
  message("NOTICE: This demo requires CGAL_Core (or LEDA), and Qt6 and will not be compiled.")
endif()
