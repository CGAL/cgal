# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Periodic_Lloyd_3_Demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()
if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

find_package(Qt5 QUIET COMPONENTS Xml Script Help OpenGL Svg)

if(Qt5Help_VERSION VERSION_LESS 5.12)
  set(CGAL_QCOLLECTIONGENERATOR_TARGET Qt5::qcollectiongenerator)
else()
  set(CGAL_QCOLLECTIONGENERATOR_TARGET Qt5::qhelpgenerator)
endif()

if(CGAL_Qt5_FOUND
   AND Qt5_FOUND
   AND TARGET ${CGAL_QCOLLECTIONGENERATOR_TARGET})

  include_directories(BEFORE ./)

  # ui file, created wih Qt Designer
  qt5_wrap_ui(uis MainWindow.ui)

  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources(CGAL_Qt5_RESOURCE_FILES ./Periodic_Lloyd_3.qrc)

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

  add_executable(
    Periodic_Lloyd_3
    Periodic_Lloyd_3.qhc
    Periodic_Lloyd_3.cpp
    MainWindow.cpp
    Viewer.cpp
    ${uis}
    ${CGAL_Qt5_RESOURCE_FILES}
    ${CGAL_Qt5_MOC_FILES})

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Periodic_Lloyd_3)

  target_link_libraries(Periodic_Lloyd_3 PRIVATE CGAL::CGAL CGAL::CGAL_Qt5
                                                 Qt5::OpenGL)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Periodic_Lloyd_3)
else(CGAL_Qt5_FOUND
  AND Qt5_FOUND
  AND QT_QCOLLECTIONGENERATOR_EXECUTABLE)

  set(PERIODIC_LLOYD_MISSING_DEPS "")

  if(NOT CGAL_Qt5_FOUND)
    set(PERIODIC_LLOYD_MISSING_DEPS
        "the CGAL Qt5 library, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()

  if(NOT Qt5_FOUND)
    set(PERIODIC_LLOYD_MISSING_DEPS "Qt5, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()

  if(NOT QT_QCOLLECTIONGENERATOR_EXECUTABLE)
    set(PERIODIC_LLOYD_MISSING_DEPS
        "${CGAL_QCOLLECTIONGENERATOR_TARGET}, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()

  message(
    STATUS
      "NOTICE: This demo requires ${PERIODIC_LLOYD_MISSING_DEPS}and will not be compiled."
  )

endif()
