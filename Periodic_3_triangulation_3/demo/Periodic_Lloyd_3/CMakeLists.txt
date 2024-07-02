# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Periodic_Lloyd_3_Demo)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

find_package(Qt6 QUIET COMPONENTS Widgets OpenGL Help ToolsTools)

set(CGAL_QCOLLECTIONGENERATOR_TARGET Qt6::qhelpgenerator)

if(CGAL_Qt6_FOUND
   AND Qt6_FOUND
   AND TARGET ${CGAL_QCOLLECTIONGENERATOR_TARGET})

  include_directories(BEFORE ./)

  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTORCC ON)

  if(DEFINED QT_QCOLLECTIONGENERATOR_EXECUTABLE)

  else()
    set(QT_QCOLLECTIONGENERATOR_EXECUTABLE qcollectiongenerator)
  endif()

  # generate QtAssistant collection file
  add_custom_command(
    OUTPUT Periodic_Lloyd_3.qhc
    DEPENDS Periodic_Lloyd_3.qhp Periodic_Lloyd_3.qhcp
    COMMAND
      ${CGAL_QCOLLECTIONGENERATOR_TARGET}
      ${CMAKE_CURRENT_SOURCE_DIR}/Periodic_Lloyd_3.qhcp -o
      ${CMAKE_CURRENT_BINARY_DIR}/Periodic_Lloyd_3.qhc
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  qt_add_executable(
    Periodic_Lloyd_3
    Periodic_Lloyd_3.qhc
    Periodic_Lloyd_3.cpp
    MainWindow.cpp
    Viewer.cpp
    Periodic_Lloyd_3.qrc)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Periodic_Lloyd_3)

  target_link_libraries(Periodic_Lloyd_3 PRIVATE CGAL::CGAL CGAL::CGAL_Qt6 Qt6::OpenGL)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Periodic_Lloyd_3)
else(CGAL_Qt6_FOUND
  AND Qt6_FOUND
  AND QT_QCOLLECTIONGENERATOR_EXECUTABLE)

  set(PERIODIC_LLOYD_MISSING_DEPS "")

  if(NOT CGAL_Qt6_FOUND)
    set(PERIODIC_LLOYD_MISSING_DEPS
        "the CGAL Qt6 library, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()

  if(NOT Qt6_FOUND)
    set(PERIODIC_LLOYD_MISSING_DEPS "Qt6, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()

  if(NOT QT_QCOLLECTIONGENERATOR_EXECUTABLE)
    set(PERIODIC_LLOYD_MISSING_DEPS
        "${CGAL_QCOLLECTIONGENERATOR_TARGET}, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()

  message("NOTICE: This demo requires ${PERIODIC_LLOYD_MISSING_DEPS} and will not be compiled.")

endif()
