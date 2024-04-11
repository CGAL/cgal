cmake_minimum_required(VERSION 3.1...3.23)
project(Periodic_4_hyperbolic_triangulation_2_Demo)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${CMAKE_BINARY_DIR})

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

find_package(CGAL REQUIRED COMPONENTS Core Qt6)

find_package(LEDA QUIET)

find_package(Qt6 QUIET COMPONENTS Widgets)
if((CGAL_Core_FOUND OR LEDA_FOUND)
   AND Qt6_FOUND
   AND CGAL_Qt6_FOUND)

  include_directories(BEFORE include)

   set(CMAKE_AUTOUIC ON)
   set(CMAKE_AUTORCC ON)

  # cpp files
  qt_add_executable(P4HDT2 P4HDT2.cpp P4HDT2.ui Main_resources.qrc)

  target_link_libraries(P4HDT2 PRIVATE CGAL::CGAL CGAL::CGAL_Qt6 Qt6::Widgets)
  if(TARGET CGAL::CGAL_Core)
    target_link_libraries(P4HDT2 PRIVATE CGAL::CGAL_Core)
  endif()
  if(LEDA_FOUND)
    target_link_libraries(P4HDT2 PRIVATE ${LEDA_LIBRARIES})
  endif()
  add_to_cached_list(CGAL_EXECUTABLE_TARGETS P4HDT2)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(P4HDT2)
else()
  message("NOTICE: This demo requires Qt6 and will not be compiled.")
endif()
