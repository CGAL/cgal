# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(GraphicsView_Demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

find_package(Qt6 QUIET COMPONENTS Widgets OpenGLWidgets Svg)

if(CGAL_Qt6_FOUND AND Qt6_FOUND)

  add_definitions(-DQT_NO_KEYWORDS)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)

  add_executable(min min.cpp ${CGAL_Qt6_RESOURCE_FILES} ${CGAL_Qt6_MOC_FILES})

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS min)

  target_link_libraries(min PRIVATE CGAL::CGAL CGAL::CGAL_Qt6 Qt6::Gui)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(min)
else()

  message("NOTICE: This demo requires CGAL and Qt6, and will not be compiled.")

endif()
