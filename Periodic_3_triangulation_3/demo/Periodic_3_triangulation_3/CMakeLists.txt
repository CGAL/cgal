# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Periodic_3_triangulation_3_Demo)

# Find CGAL
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

# Find Qt6 itself
find_package(Qt6 QUIET COMPONENTS Widgets OpenGL Help ToolsTools)

if(Qt6_FOUND)
  add_definitions(-DQT_NO_KEYWORDS)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)
endif(Qt6_FOUND)

set(CGAL_QCOLLECTIONGENERATOR_TARGET Qt6::qhelpgenerator)

if(CGAL_Qt6_FOUND
   AND Qt6_FOUND
   AND TARGET ${CGAL_QCOLLECTIONGENERATOR_TARGET})

   set(CMAKE_AUTOMOC ON)
   set(CMAKE_AUTOUIC ON)
   set(CMAKE_AUTORCC ON)

  # generate QtAssistant collection file
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Periodic_3_triangulation_3.qhc
    DEPENDS Periodic_3_triangulation_3.qhp Periodic_3_triangulation_3.qhcp
    COMMAND
      ${CGAL_QCOLLECTIONGENERATOR_TARGET}
      ${CMAKE_CURRENT_SOURCE_DIR}/Periodic_3_triangulation_3.qhcp -o
      ${CMAKE_CURRENT_BINARY_DIR}/Periodic_3_triangulation_3.qhc)

  # The executable itself
  qt_add_executable(
    periodic_3_triangulation_3_demo
    Scene.cpp
    Viewer.cpp
    MainWindow.cpp
    periodic_3_triangulation_3_demo.cpp
    MainWindow.ui
    Periodic_3_triangulation_3.qhc
    Periodic_3_triangulation_3.qrc)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS periodic_3_triangulation_3_demo)

  target_link_libraries(periodic_3_triangulation_3_demo
                        PRIVATE CGAL::CGAL CGAL::CGAL_Qt6 Qt6::OpenGL)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(periodic_3_triangulation_3_demo)
else()

  set(PERIODIC_TRIANGULATION_MISSING_DEPS "")

  if(NOT CGAL_Qt6_FOUND)
    set(PERIODIC_TRIANGULATION_MISSING_DEPS
        "the CGAL Qt6 library, ${PERIODIC_TRIANGULATION_MISSING_DEPS}")
  endif()

  if(NOT Qt6_FOUND)
    set(PERIODIC_TRIANGULATION_MISSING_DEPS
        "Qt6, ${PERIODIC_TRIANGULATION_MISSING_DEPS}")
  endif()

  if(NOT TARGET ${CGAL_QCOLLECTIONGENERATOR_TARGET})
    set(PERIODIC_TRIANGULATION_MISSING_DEPS
        "${CGAL_QCOLLECTIONGENERATOR_TARGET}, ${PERIODIC_TRIANGULATION_MISSING_DEPS}"
    )
  endif()

  message("NOTICE: This demo requires ${PERIODIC_TRIANGULATION_MISSING_DEPS}and will not be compiled.")

endif()
