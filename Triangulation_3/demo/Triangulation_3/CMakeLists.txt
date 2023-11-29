# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Triangulation_3_Demo)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

find_package(Qt6 QUIET COMPONENTS Widgets OpenGL)

if(Qt6_FOUND)
  add_definitions(-DQT_NO_KEYWORDS)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)
endif(Qt6_FOUND)

# Activate concurrency ? (turned OFF by default)
option(CGAL_ACTIVATE_CONCURRENT_TRIANGULATION_3
       "Activate parallelism in Triangulation_3" OFF)

# And add -DCGAL_CONCURRENT_TRIANGULATION_3 if that option is ON
if(CGAL_ACTIVATE_CONCURRENT_TRIANGULATION_3)
  add_definitions(-DCGAL_CONCURRENT_TRIANGULATION_3)
  find_package(TBB REQUIRED)
  include(CGAL_TBB_support)
else(CGAL_ACTIVATE_CONCURRENT_TRIANGULATION_3)
  option(LINK_WITH_TBB
         "Link with TBB anyway so we can use TBB timers for profiling" ON)
  if(LINK_WITH_TBB)
    find_package(TBB)
    include(CGAL_TBB_support)
  endif(LINK_WITH_TBB)
endif()

if(CGAL_Qt6_FOUND AND Qt6_FOUND)

  include_directories(BEFORE ./)
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTORCC ON)

  # cpp files
  qt_add_executable(
    T3_demo
    T3_demo.cpp
    MainWindow.cpp
    Viewer.cpp
    PreferenceDlg.cpp
    Scene.cpp
    MainWindow.ui
    T3_demo.qrc)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS T3_demo)

  target_link_libraries(T3_demo PRIVATE CGAL::CGAL CGAL::CGAL_Qt6)
  target_link_libraries(T3_demo PRIVATE Qt6::OpenGL)
  if(TARGET CGAL::TBB_support)
    target_link_libraries(T3_demo PUBLIC CGAL::TBB_support)
  endif()

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(T3_demo)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(T3_demo)

else(CGAL_Qt6_FOUND AND Qt6_FOUND)

  set(TRIANGULATION_3_MISSING_DEPS "")

  if(NOT CGAL_Qt6_FOUND)
    set(TRIANGULATION_3_MISSING_DEPS
        "the CGAL Qt6 library, ${TRIANGULATION_3_MISSING_DEPS}")
  endif()

  if(NOT Qt6_FOUND)
    set(TRIANGULATION_3_MISSING_DEPS "Qt6, ${TRIANGULATION_3_MISSING_DEPS}")
  endif()

  message("NOTICE: This demo requires ${TRIANGULATION_3_MISSING_DEPS}, and will not be compiled.")

endif(CGAL_Qt6_FOUND AND Qt6_FOUND)
