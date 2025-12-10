# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Alpha_wrap_2_Demo)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

find_package(Qt6 QUIET COMPONENTS Widgets)

if(CGAL_Qt6_FOUND AND Qt6_FOUND)
  add_compile_definitions(QT_NO_KEYWORDS)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTORCC ON)

  # The executable itself.
  qt6_wrap_ui(Alpha_wrap_2_uis Alpha_wrap_2.ui Alpha_wrap_2_options.ui Screenshot_options.ui)
  qt6_wrap_cpp(Alpha_wrap_2_mocs Alpha_wrap_2_options.h Screenshot_options.h)
  qt6_add_resources(Alpha_wrap_2_qrc Alpha_wrap_2.qrc)

  qt_add_executable(Alpha_wrap_2 Alpha_wrap_2.cpp Viewer.cpp ${Alpha_wrap_2_uis} ${Alpha_wrap_2_mocs} ${Alpha_wrap_2_qrc})

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Alpha_wrap_2)

  target_link_libraries(Alpha_wrap_2 PRIVATE CGAL::CGAL CGAL::CGAL_Qt6 Qt6::Widgets)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Alpha_wrap_2)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Alpha_wrap_2)
else()
  message("NOTICE: This demo requires CGAL and Qt6, and will not be compiled.")
endif()
